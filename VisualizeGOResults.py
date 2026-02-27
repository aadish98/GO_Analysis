#!/usr/bin/env python3
import argparse, pathlib, re, math
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import to_rgb, LinearSegmentedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colorbar import ColorbarBase
import matplotlib.patheffects as pe
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch

# ────────────────────────── tweakables ──────────────────────────
ROOT = pathlib.Path(__file__).resolve().parent
CAT_FILE = ROOT / "Data" / "GO-Term_categories.xlsx"
CLUSTER_FILE = ROOT / "Data" / "CSW_7datasets" / "CSW_7datasets_AllTerms_Mapped.xlsx"

TOP_N, MAX_ROWS = 1000, 100000
MAX_TERMS_PER_PLOT = 15  # Maximum terms per plot before splitting
X_PAD = 2.0
LABEL_FACTOR = 1.0
EXTRA_PAD = 1
TEXT_OUTLINE = [pe.withStroke(linewidth=3, foreground="white")]
MIN_LEGEND_HEIGHT = 6.0
# ────────────────────────────────────────────────────────────────

# ─── keyword & colour table ─────────────────────────────────────
cat_df = pd.read_excel(CAT_FILE)
cat_df["Keywords"] = cat_df["Keywords"].str.split(",").apply(
    lambda L: [k.strip().lower() for k in L])
GROUPS = cat_df["Group"].tolist()
GROUP2KW   = dict(zip(GROUPS, cat_df["Keywords"]))
GROUP2HEX  = dict(zip(GROUPS, cat_df["Hexcode"]))
GROUP2LABEL= {g: g for g in GROUPS}
ASPECT_HEX = {
    "Biological Process":   "#1f77b4",  # blue
    "Molecular Function":   "#ff7f0e",  # orange
    "Cellular Component":   "#2ca02c",  # green
    "Kegg Pathway": "#9467bd",  # purple
}
NO_BROAD_GROUP_PREFIX = "__NO_BROAD__::"


def make_no_broad_group(term):
    term_label = pretty(str(term)).strip() if term is not None else ""
    if not term_label:
        term_label = "Unknown Term"
    return f"{NO_BROAD_GROUP_PREFIX}{term_label}"


def is_no_broad_group(group):
    return str(group).startswith(NO_BROAD_GROUP_PREFIX)

def aspect_text_color(aspect):
    return ASPECT_HEX.get(aspect, "#333333")

def add_content_box(ax, artists, *, pad_x=0.02, pad_y=0.012, min_w=0.30, min_h=0.05, zorder=0):
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bboxes = []
    for artist in artists:
        if artist is None:
            continue
        try:
            bb = artist.get_window_extent(renderer=renderer)
        except Exception:
            continue
        if bb is None:
            continue
        vals = np.array([bb.x0, bb.y0, bb.x1, bb.y1], dtype=float)
        if not np.all(np.isfinite(vals)):
            continue
        bboxes.append(bb)

    if not bboxes:
        return None

    x0 = min(bb.x0 for bb in bboxes)
    y0 = min(bb.y0 for bb in bboxes)
    x1 = max(bb.x1 for bb in bboxes)
    y1 = max(bb.y1 for bb in bboxes)

    inv = ax.transAxes.inverted()
    ax_x0, ax_y0 = inv.transform((x0, y0))
    ax_x1, ax_y1 = inv.transform((x1, y1))

    left = max(0.01, min(ax_x0, ax_x1) - pad_x)
    bottom = max(0.01, min(ax_y0, ax_y1) - pad_y)
    right = min(0.99, max(ax_x0, ax_x1) + pad_x)
    top = min(0.99, max(ax_y0, ax_y1) + pad_y)

    width = max(min_w, right - left)
    height = max(min_h, top - bottom)
    width = min(width, 0.99 - left)
    height = min(height, 0.99 - bottom)

    patch = FancyBboxPatch(
        (left, bottom), width, height,
        boxstyle="round,pad=0.01", transform=ax.transAxes,
        facecolor="white", edgecolor="black", linewidth=0.9, zorder=zorder
    )
    ax.add_patch(patch)
    return patch

assign_group = lambda term: next((g for g,kws in GROUP2KW.items()
                                  if any(k in term.lower() for k in kws)),None)
pretty = lambda s: re.sub(r'^.*?:', '', s).split('~')[-1].strip()
darken = lambda h,f: tuple(np.array(to_rgb(h))*f)

# ─── Cluster Assignment & Hierarchical Packing ──────────────────
class ClusterAssigner:
    def __init__(self, path, unknown_policy="jaccard", source_name=None):
        """
        unknown_policy:
          - "jaccard": assign unknown terms to best existing group by Jaccard (based on Genes)
          - "unique": create a new numeric group for every unknown term (one group per term)
        """
        self.unknown_policy = unknown_policy
        self._unknown_term_to_group = {}  # normalized_term_key -> int group id
        self.source_name = source_name
        self.source_sheet = None
        self.no_broad_term_keys = set()

        # ---------- helpers ----------
        # Use your pretty() function to normalize term names so matching is consistent
        self._term_key = lambda t: pretty(str(t)).strip().lower()
        self._name_key = lambda s: re.sub(r"[^a-z0-9]+", "", pathlib.Path(str(s)).stem.lower())

        # ---------- load cluster file (CSV or Excel) ----------
        if not path.exists():
            print(f"Warning: {path} not found. Clustering disabled.")
            self.df = pd.DataFrame(columns=["Group", "Term", "Genes"])
        else:
            try:
                if path.suffix.lower() == ".csv":
                    self.df = pd.read_csv(path)
                else:
                    raw = pd.read_excel(path, sheet_name=None)
                    if isinstance(raw, dict):
                        if source_name:
                            source_key = self._name_key(source_name)
                            sheet_match = next(
                                (sheet for sheet in raw if self._name_key(sheet) == source_key),
                                None
                            )
                            if sheet_match is not None:
                                self.source_sheet = sheet_match
                                self.df = raw[sheet_match].copy()
                            else:
                                available = ", ".join(map(str, raw.keys()))
                                print(
                                    f"Warning: no matching sheet in {path.name} for source '{source_name}'. "
                                    f"Falling back to all sheets. Available sheets: {available}"
                                )
                                self.df = pd.concat(raw.values(), ignore_index=True)
                        else:
                            self.df = pd.concat(raw.values(), ignore_index=True)
                    else:
                        self.df = raw

                # normalize col names — support both canonical and CSV-style headers
                cols = {c.lower().strip().replace("_", " "): c for c in self.df.columns}
                rename = {}
                if "term" in cols:       rename[cols["term"]]       = "Term"
                if "go term" in cols:    rename[cols["go term"]]    = "Term"
                if "group" in cols:      rename[cols["group"]]      = "Group"
                if "broad category" in cols: rename[cols["broad category"]] = "Group"
                if "genes" in cols:      rename[cols["genes"]]      = "Genes"
                self.df = self.df.rename(columns=rename)

                # If multiple source columns map to the same canonical name (e.g., Term + GO_Term),
                # coalesce them into one Series so downstream .str operations stay valid.
                for canonical in ["Term", "Group", "Genes"]:
                    dup_idx = [i for i, c in enumerate(self.df.columns) if c == canonical]
                    if len(dup_idx) > 1:
                        merged = self.df.iloc[:, dup_idx[0]]
                        for j in dup_idx[1:]:
                            merged = merged.combine_first(self.df.iloc[:, j])
                        keep_pos = [True] * len(self.df.columns)
                        for j in dup_idx[1:]:
                            keep_pos[j] = False
                        self.df = self.df.iloc[:, keep_pos].copy()
                        self.df[canonical] = merged

                for c in ["Term", "Group"]:
                    if c not in self.df.columns:
                        raise ValueError(f"{path.name} must contain a column mappable to '{c}'.")

                group_raw = self.df["Group"]
                self.df["Term"] = self.df["Term"].astype(str).str.strip()
                missing_group_mask = group_raw.isna() | group_raw.astype(str).str.strip().eq("")
                self.no_broad_term_keys = set(
                    self.df.loc[missing_group_mask, "Term"]
                    .astype(str)
                    .map(self._term_key)
                    .tolist()
                )

                self.df["Group"] = group_raw.astype(str).str.strip()
                self.df = self.df[(self.df["Group"] != "") & (self.df["Term"] != "")].copy()

                if "Genes" not in self.df.columns:
                    self.df["Genes"] = ""

            except Exception as e:
                print(f"Error loading {path}: {e}")
                self.df = pd.DataFrame(columns=["Group", "Term", "Genes"])
                self.no_broad_term_keys = set()

        if self.source_sheet is not None:
            self.assignment_scope = f"{path.name} [sheet: {self.source_sheet}]"
        else:
            self.assignment_scope = path.name

        # ---------- build term->group mapping (robust) ----------
        # IMPORTANT: normalize term keys so terms defined in xlsx don't get treated as "unknown"
        self.term_to_group = {}
        for _, row in self.df.iterrows():
            k = self._term_key(row["Term"])

            # cast group to int if possible, else keep as stripped string
            g_raw = row["Group"]
            g_num = pd.to_numeric(g_raw, errors="coerce")
            g_val = int(g_num) if pd.notna(g_num) else str(g_raw).strip()

            # first occurrence wins (so you can override by ordering sheets if desired)
            if k not in self.term_to_group:
                self.term_to_group[k] = g_val

        # ---------- for jaccard mode: group -> union genes ----------
        self.group_genes = {}
        for _, row in self.df.iterrows():
            g_raw = row["Group"]
            g_num = pd.to_numeric(g_raw, errors="coerce")
            g_val = int(g_num) if pd.notna(g_num) else str(g_raw).strip()

            gs = self._parse_genes(row.get("Genes", ""))
            self.group_genes.setdefault(g_val, set()).update(gs)

        # ---------- determine next integer for NEW groups ----------
        # NEW numeric ids should start after max numeric group already present
        existing_numeric = pd.to_numeric(self.df["Group"], errors="coerce")
        if existing_numeric.notna().any():
            self._next_new_id = int(existing_numeric.max()) + 1
        else:
            self._next_new_id = 1

        # ---------- known groups list (for your GROUPS init) ----------
        # keep stable ordering: numeric ascending, then strings
        known = list(self.group_genes.keys())
        self.known_groups = sorted(known, key=lambda x: (isinstance(x, str), x))

    def _parse_genes(self, x):
        if pd.isna(x):
            return set()
        return set(re.split(r'[,;\s\n]+', str(x).strip())) - {""}

    def _make_new_group_for_term(self, term_key):
        # stable mapping within a run
        if term_key in self._unknown_term_to_group:
            return self._unknown_term_to_group[term_key]
        gid = int(self._next_new_id)
        self._next_new_id += 1
        self._unknown_term_to_group[term_key] = gid
        return gid

    def normalize_term_key(self, term):
        return self._term_key(term)

    def get_group(self, term, genes_val):
        # NaN/empty term fallback
        if pd.isna(term) or str(term).strip() == "":
            return "Unassigned"

        k = self._term_key(term)

        # Known term from xlsx (robust key)
        if k in self.term_to_group:
            return self.term_to_group[k]

        # Unique mode: one new integer group per unknown term
        if self.unknown_policy == "unique":
            return self._make_new_group_for_term(k)

        # Jaccard mode: assign to best existing group by gene overlap
        t_genes = self._parse_genes(genes_val)
        if not t_genes:
            return "Unassigned"

        best_g = "Unassigned"
        best_j = -1.0
        for g, g_genes in self.group_genes.items():
            if not g_genes:
                continue
            inter = len(t_genes & g_genes)
            union = len(t_genes | g_genes)
            j = inter / union if union > 0 else 0.0
            if j > best_j:
                best_j = j
                best_g = g
        return best_g

def get_jaccard(genes_a, genes_b):
    sa = set(re.split(r'[,;\s\n]+', str(genes_a).strip())) - {''}
    sb = set(re.split(r'[,;\s\n]+', str(genes_b).strip())) - {''}
    if not sa or not sb: return 0.0
    return len(sa & sb) / len(sa | sb)

def cluster_terms_greedy(df_group):
    # Sort by significance (most significant first -> smallest Pvalue)
    df_sorted = df_group.sort_values("Pvalue", ascending=True)
    
    clusters = [] # List of lists of indices
    processed = set()
    
    for idx, row in df_sorted.iterrows():
        if idx in processed:
            continue
            
        current_cluster = [idx]
        processed.add(idx)
        genes_a = row.get("Genes", "")
        
        # Look for other terms that match
        for other_idx, other_row in df_sorted.iterrows():
            if other_idx in processed:
                continue
            
            genes_b = other_row.get("Genes", "")
            if get_jaccard(genes_a, genes_b) > 0.7:
                current_cluster.append(other_idx)
                processed.add(other_idx)
        
        clusters.append(current_cluster)
        
    return clusters

def pack_layout_clustered(df, step):
    # Returns:
    # label_y_map: {pretty_term: y}
    # dot_pos_map: {pretty_term: (x, y)}
    # cluster_info: list of dicts (for CSV export)
    
    label_y_map = {}
    dot_pos_map = {}
    cluster_info = []
    cur_y = 0.0
    
    # Process groups in order consistent with aspect_plot
    # We iterate GROUPS reversed to match pack_y_aspect logic (bottom-up packing)
    
    present_groups = [g for g in GROUPS if g in df["Group"].unique()]
    # Reverse them to pack bottom-up
    loop_groups = present_groups[::-1]
    
    for g in loop_groups:
        sub = df[df["Group"] == g]
        if sub.empty: continue
        
        clusters = cluster_terms_greedy(sub)
        
        # clusters are ordered Most Sig -> Least Sig.
        # We are packing Bottom -> Top.
        # We want Most Sig at Top of group.
        # So we should process Least Sig clusters first.
        clusters_reversed = clusters[::-1]
        
        for clus_indices in clusters_reversed:
            clus_rows = sub.loc[clus_indices]
            
            # Representative: Most significant term (first in clus_indices)
            rep_idx = clus_indices[0]
            rep_row = sub.loc[rep_idx]
            rep_x = rep_row["-log10P"]
            rep_term = rep_row["Term"]
            
            # Record cluster info if it's an "overlapping" dot (size > 1)
            if len(clus_indices) > 1:
                for idx in clus_indices:
                    row = sub.loc[idx]
                    cluster_info.append({
                        "Aspect": row.get("Aspect", "Unknown"),
                        "Group": g,
                        "Cluster_Representative": rep_term,
                        "Term": row["Term"],
                        "Pvalue": row["Pvalue"],
                        "Genes": row.get("Genes", ""),
                        "Cluster_Size": len(clus_indices)
                    })

            # Sort terms in cluster for label stacking
            clus_rows_sorted = clus_rows.sort_values("Pvalue", ascending=False)
            
            term_ys = []
            for _, r in clus_rows_sorted.iterrows():
                term_p = r["pretty"]
                label_y_map[term_p] = cur_y
                term_ys.append(cur_y)
                cur_y += step
            
            # Dot Y: Average of term_ys
            dot_y = sum(term_ys) / len(term_ys)
            
            for _, r in clus_rows.iterrows():
                dot_pos_map[r["pretty"]] = (rep_x, dot_y)
            
        cur_y += step * 0.6 # Gap between groups
        
    return label_y_map, dot_pos_map, cluster_info

def pack_y_hierarchical(df, step):
    # This function is deprecated in favor of aspect_plot + pack_y_aspect, 
    # but kept if needed for backward compatibility or single-panel mode
    return pack_y_aspect(df, step)

def pack_y_aspect(df, step):
    y_map = {}
    cur_y = 0.0
    # Sort groups reversed so first processed gets lowest Y? 
    # We want top of list to be top of plot. 
    # Usually scatter plots have Y=0 at bottom.
    # So first item should be at top (highest Y).
    # If we iterate top-down, we decrease Y.
    # Or iterate bottom-up, we increase Y.
    # Let's iterate Groups in standard order (0, 1, 2...)
    # But within group, sort by P-value ascending (most sig first? no, -logP is X)
    # The requirement: "sorted by Group assignment... and then sorted by -log10P value (descending)"
    # Descending significance = largest -log10P first.
    # Usually top item in list is most significant.
    # So we want most significant at Y_max.
    # So we should process:
    #   Group 0 (top)
    #      Term 1 (most sig) -> Y_max
    #      Term 2            -> Y_max - 1
    #   Group 1
    # ...
    
    # To assign Y coordinates starting from 0 (bottom) going up:
    # We should process groups in REVERSE order (bottom group first).
    # Within group, process LEAST significant first (bottom term first).
    
    groups = sorted(df["Group"].dropna().unique(), reverse=True) # Bottom group first?
    # Actually, if we want Group 0 at top, and we are packing from 0 upwards:
    # We should process Group N first (at Y=0), then Group N-1... up to Group 0 (at Y=Max).
    # Yes.
    
    for group in groups:
        # We want most significant at top of group.
        # So we process least significant first (bottom of group).
        # Sort by Pvalue descending (largest P = smallest -logP = least significant)
        grp_df = df[df["Group"] == group].sort_values("Pvalue", ascending=False)
        
        last_genes = None
        last_y = cur_y
        
        for _, row in grp_df.iterrows():
            overlap = False
            if last_genes is not None:
                sim = get_jaccard(row.get("Genes"), last_genes)
                if sim > 0.7:
                    overlap = True
            
            if overlap:
                y_map[row["pretty"]] = last_y
                last_genes = row.get("Genes")
            else:
                cur_y += step
                y_map[row["pretty"]] = cur_y
                last_y = cur_y
                last_genes = row.get("Genes")
        
        cur_y += step * 0.6 # Gap between groups
        
    return y_map


def pack_y_grouped(df, step):
    y, cur = {}, 0.0
    for g in GROUPS:
        sub = df[df["Group"]==g]
        for _,row in sub.iterrows():
            y[row["pretty"]] = cur; cur += step
        if not sub.empty: cur += step*.6
    return y

# ─── XLSX loader ────────────────────────────────────────────────
def read_csv(p, *, keep_uncat=False, aspect_mode=False, cluster_assigner=None):
    # Load Excel; if multiple sheets exist, concatenate them
    raw = pd.read_excel(p, sheet_name=None)
    if isinstance(raw, dict):
        df = pd.concat(raw.values(), ignore_index=True)
    else:
        df = raw
    df = df.rename(columns=str.strip)
    df["%"]      = pd.to_numeric(df["%"], errors="coerce")
    df["Pvalue"] = pd.to_numeric(df["Pvalue"], errors="coerce")
    df["-log10P"] = -np.log10(df["Pvalue"].clip(lower=1e-300))

    if not aspect_mode:
        df["Group"] = df["Term"].apply(assign_group)
    else:
        # ── map each row to high-level aspect when --by-aspect is active
        def _aspect(cat):  # cat = value from the Category column
            if pd.isna(cat):
                return None
            s = str(cat).upper()
            if "GOTERM_BP" in s or "BIOLOGICAL_PROCESS" in s:
                return "Biological Process"
            if "GOTERM_CC" in s or "CELLULAR_COMPONENT" in s:
                return "Cellular Component"
            if "GOTERM_MF" in s or "MOLECULAR_FUNCTION" in s:
                return "Molecular Function"
            if "KEGG" in s or "PATHWAY" in s:
                return "Kegg Pathway"
            return None

        # Category must exist – raise a helpful error if not
        if "Category" not in df.columns:
            raise ValueError(f"{p.name} lacks a 'Category' column required for --by_GO")

        df["Aspect"] = df["Category"].map(_aspect)
        
        # Use ClusterAssigner for Group
        if cluster_assigner:
            term_keys = df["Term"].map(cluster_assigner.normalize_term_key)
            known_keys = set(cluster_assigner.term_to_group.keys())
            no_broad_mask = term_keys.isin(cluster_assigner.no_broad_term_keys)
            missing_mask = ~term_keys.isin(known_keys) & ~no_broad_mask
            if missing_mask.any():
                unmapped_terms = (df.loc[missing_mask, "Term"]
                                  .map(pretty)
                                  .dropna()
                                  .drop_duplicates()
                                  .tolist())
                preview_n = 25
                shown = unmapped_terms[:preview_n]
                more = len(unmapped_terms) - len(shown)
                listing = "\n".join(f"  - {t}" for t in shown)
                if more > 0:
                    listing += f"\n  ... and {more} more unmapped terms."
                print(
                    "\n"
                    "Warning: unmapped GO terms detected; assigning them to unique fallback groups.\n"
                    f"Mapping source: {cluster_assigner.assignment_scope}\n"
                    f"Unmapped term count: {len(unmapped_terms)}\n"
                    f"{listing}\n"
                )
            df["NoBroadCategory"] = no_broad_mask
            df["Group"] = df.apply(
                lambda row: make_no_broad_group(row.get("Term"))
                if bool(row.get("NoBroadCategory", False))
                else cluster_assigner.get_group(row.get("Term"), row.get("Genes")),
                axis=1,
            )
        else:
             df["Group"] = "Unassigned"

    if not keep_uncat:
        if aspect_mode:
             df = df.dropna(subset=["Aspect"])
        else:
             df = df.dropna(subset=["Group"])
    else:
        if aspect_mode:
            df["Aspect"] = df["Aspect"].fillna("Other")
            if "Unassigned" not in GROUPS:
                 # Logic for Unassigned group handled in main/bubble_plot
                 pass
        else:
            df["Group"] = df["Group"].fillna("Other")
            if "Other" not in GROUPS:
                GROUPS.append("Other")
                GROUP2HEX["Other"]   = "#7f7f7f"
                GROUP2KW ["Other"]   = []
                GROUP2LABEL["Other"] = "Other"

    return (df.sort_values("Pvalue")
              .drop_duplicates("Term", keep="first"))

def validate_all_terms_mapped(xlsx_paths, cluster_assigner):
    known_keys = set(cluster_assigner.term_to_group.keys())
    missing_by_file = {}
    for xlsx_path in xlsx_paths:
        raw = pd.read_excel(xlsx_path, sheet_name=None)
        if isinstance(raw, dict):
            df = pd.concat(raw.values(), ignore_index=True)
        else:
            df = raw
        df = df.rename(columns=str.strip)
        if "Term" not in df.columns:
            continue
        keys = df["Term"].astype(str).map(cluster_assigner.normalize_term_key)
        no_broad_mask = keys.isin(cluster_assigner.no_broad_term_keys)
        missing_mask = ~keys.isin(known_keys) & ~no_broad_mask
        if missing_mask.any():
            missing_terms = (df.loc[missing_mask, "Term"]
                             .map(pretty)
                             .dropna()
                             .drop_duplicates()
                             .tolist())
            if missing_terms:
                missing_by_file[xlsx_path.name] = missing_terms

    if missing_by_file:
        total = sum(len(v) for v in missing_by_file.values())
        lines = []
        for file_name, terms in missing_by_file.items():
            sample = ", ".join(terms[:8])
            more = len(terms) - min(8, len(terms))
            suffix = f", ... (+{more} more)" if more > 0 else ""
            lines.append(f"  - {file_name}: {sample}{suffix}")
        details = "\n".join(lines)
        print(
            "\n"
            "Warning: GO terms missing from mapping file. "
            "These terms will be assigned fallback unique groups.\n"
            f"Mapping source: {cluster_assigner.assignment_scope}\n"
            f"Unmapped term count: {total}\n"
            f"Files with unmapped terms:\n{details}\n"
        )
    return missing_by_file

# ─── plotting ───────────────────────────────────────────────────
def bubble_plot(df, png=None, *, split=False, side=None,
                fontsize=16, spacing=10000, skip_set="", ax=None, plot_width=3):
    # skip terms?
    if skip_set:
        df = df[~df["Term"].map(pretty).str.lower().isin(skip_set)]

    # reflect wake
    if split and side == "wake":
        df = df.copy();
        df["-log10P"] *= -1

    # top-N per group
    keep = [df[df["Group"] == g].nsmallest(TOP_N, "Pvalue") for g in GROUPS]
    top = pd.concat(keep, ignore_index=True)
    if top.empty: return
    top["pretty"] = top["Term"].map(pretty)
    top = (top.sort_values(["Group", "Pvalue"])
           .drop_duplicates("pretty", keep="first")
           .head(MAX_ROWS))
    top["label"] = top["pretty"]

    # enforce panel side
    if side is not None:
        delta = .4
        if side == "sleep":
            top.loc[top["-log10P"] < 0, "-log10P"] *= -1
        else:
            top.loc[top["-log10P"] > 0, "-log10P"] *= -1
        top["-log10P"] += delta * np.sign(top["-log10P"])

    # vertical packing
    row_step = max(.8, fontsize / 12 * .35)
    if 'side' in top.columns:
        y_map = {**pack_y_grouped(top[top.side == 'wake'], row_step),
                 **pack_y_grouped(top[top.side == 'sleep'], row_step)}
    elif 'Aspect' in top.columns:
        y_map = pack_y_hierarchical(top, row_step)
    else:
        y_map = pack_y_grouped(top, row_step)
    if not y_map: return
    top["y"] = top["pretty"].map(y_map)
    max_y = max(y_map.values())

    # horizontal limits
    # Include the extra padding in the limit calculation
    extra_pad = 1.0
    half_span = max(abs(p) + LABEL_FACTOR * len(lbl) + (X_PAD + extra_pad)
                    for p, lbl in zip(top["-log10P"], top["label"]))
    xlim_spacer = 15 if spacing >= 10000 else 15
    xlim = int(math.ceil((half_span) / 10.0) * xlim_spacer)

    fig_w = max(12, xlim / plot_width);
    fig_h = max(MIN_LEGEND_HEIGHT, 0.45 * (max_y + 2.5))

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 0.25], wspace=0.1)
    ax = fig.add_subplot(gs[0, 0])
    legend_ax = fig.add_subplot(gs[0, 1])

    # axes range
    if split:
        ax.set_xlim(-xlim, xlim);
        ax.axvline(0, color='steelblue', lw=2, zorder=2)
    else:
        ax.set_xlim(0, xlim);
        ax.axvline(0, color='steelblue', lw=2, zorder=2)

    scale = 650 / np.sqrt(top["%"].max() or 1)

    # ── draw bubbles & labels + group brackets ──────────────────
    for g in GROUPS:
        sub = top[top["Group"] == g]
        if sub.empty: continue
        base = GROUP2HEX[g]
        
        ax.scatter(sub["-log10P"], sub["y"],
                   s=np.sqrt(sub["%"]) * scale, color=base,
                   edgecolor='black', lw=.6, zorder=3)

        for x, y, lbl, aspect in zip(sub["-log10P"], sub["y"], sub["label"], sub.get("Aspect", pd.Series([""] * len(sub)))):
            # Increased padding for normal plots (from X_PAD to X_PAD + 1)
            extra_pad = 1.0 
            dx, ha = (-(X_PAD + extra_pad), 'right') if (split and x < 0) else ((X_PAD + extra_pad), 'left')
            ax.text(x + dx, y, lbl, ha=ha, va='center', fontsize=fontsize * .8,
                    color=aspect_text_color(aspect), zorder=4, path_effects=TEXT_OUTLINE)

        # ─── bracket & label
        sub_brackets = []
        if 'Aspect' in top.columns:
            for asp, asp_df in sub.groupby("Aspect"):
                sub_brackets.append(asp_df)
        else:
            sub_brackets.append(sub)

        for s_b in sub_brackets:
            y0, y1 = s_b["y"].min() - 0.4, s_b["y"].max() + 0.4
            mid = (y0 + y1) / 2
            pad_x = 0.2

            if split and s_b["-log10P"].mean() < 0:  # --- WAKE (left side)
                # Include extra padding in bracket calculation
                extra_pad = 1.0
                text_edges = (s_b["-log10P"] - (X_PAD + extra_pad)
                              - LABEL_FACTOR * s_b["label"].str.len())
                xb = text_edges.min() - pad_x
                ax.plot([xb, xb], [y0, y1], color=base, lw=1.2)
                ax.plot([xb, xb + 1], [y0, y0], color=base, lw=1.2)
                ax.plot([xb, xb + 1], [y1, y1], color=base, lw=1.2)
                ax.text(xb - 0.8, mid, GROUP2LABEL[g],
                        ha='right', va='center',
                        fontsize=fontsize * .85, color=base, weight='bold')

            else:  # --- SLEEP or single-panel
                # Include extra padding in bracket calculation
                extra_pad = 1.0
                text_edges = (s_b["-log10P"] + (X_PAD + extra_pad)
                              + LABEL_FACTOR * s_b["label"].str.len())
                xb = text_edges.max() + pad_x
                ax.plot([xb, xb], [y0, y1], color=base, lw=1.2)
                ax.plot([xb - 1, xb], [y0, y0], color=base, lw=1.2)
                ax.plot([xb - 1, xb], [y1, y1], color=base, lw=1.2)
                ax.text(xb + 0.8, mid, GROUP2LABEL[g],
                        ha='left', va='center',
                        fontsize=fontsize * .85, color=base, weight='bold')

    # ─── Aspect Headers
    if 'Aspect' in top.columns:
        for aspect in top['Aspect'].unique():
            sub_a = top[top['Aspect'] == aspect]
            if sub_a.empty: continue
            ys = sub_a["y"]
            min_y, max_y = ys.min(), ys.max()
            mid_y = (min_y + max_y) / 2
            
            # Position text outside right edge
            ax.text(1.02, mid_y, aspect, ha='left', va='center', 
                    fontsize=fontsize+2, weight='bold', rotation=270, 
                    transform=ax.get_yaxis_transform(), color='#444444')

    # ────────────────────────────────────────────────────────────

    # MODIFICATION: Reduced top padding from +2 to +1
    ax.set_ylim(-1, max_y + 1);
    ax.set_yticks([])

    if split:
        ax.set_xticks(np.arange(-xlim, xlim + 1, 5))
        ax.set_xticklabels([abs(t) for t in ax.get_xticks()], fontsize=fontsize * .9)
    else:
        ax.set_xticks(np.arange(0, xlim + 1, 5))
        ax.set_xticklabels(ax.get_xticks(), fontsize=fontsize * .9)

    ax.grid(axis='x', ls='--', alpha=.15)
    ax.set_xlabel("Negative log p-value", fontsize=fontsize + 2, weight='bold')

    if split:
        hdr = max_y + 1.5
        ax.text(-xlim + 1, hdr, "Wake Genes GO Terms", ha='left', va='top',
                weight='bold', fontsize=fontsize + 2)
        ax.text(xlim - 1, hdr, "Sleep Genes GO Terms", ha='right', va='top',
                weight='bold', fontsize=fontsize + 2)

    legend_ax.axis("off")
    present = [g for g in GROUPS if g in top["Group"].values][::-1]
    if not present:
        fig.tight_layout()
        fig.savefig(png, dpi=300, bbox_inches="tight", pad_inches=0.35);
        plt.close(fig);
        return

    handles = [plt.Line2D([0], [0], marker='o', ls='', color=GROUP2HEX[g],
                          markersize=8) for g in present]
    leg = legend_ax.legend(handles, [GROUP2LABEL[g] for g in present],
                           title="Group", loc="upper left", frameon=False,
                           fontsize=fontsize * .75, title_fontsize=fontsize * .8)
    legend_ax.add_artist(leg);
    fig.canvas.draw()
    y_cur = leg.get_window_extent().transformed(legend_ax.transAxes.inverted()).y0

    t_fs, label_fs = fontsize * .8, fontsize * .75

    # MODIFICATION START: Revised legend spacing for better layout
    y_cur -= 0.10  # Space after Group list
    legend_ax.text(.5, y_cur, "% Genes", ha='center', va='bottom',
                   fontsize=t_fs, transform=legend_ax.transAxes)

    y_cur -= 0.04 # Space after % Genes title
    sizes_to_show = [100, 75, 50, 25, 10, 5]
    x_b, x_l = .35, .65
    fig.canvas.draw()
    bottom_limit = 0.02
    while True:
        s_vals, row_step, title_gap = _size_legend_layout(fig, legend_ax, sizes_to_show, scale)
        first_y = y_cur - title_gap
        last_y = first_y - (len(sizes_to_show) - 1) * row_step
        if last_y >= bottom_limit or len(sizes_to_show) <= 2:
            break
        sizes_to_show = sizes_to_show[:-1]
    y_cur = first_y
    for pct, s_val in zip(sizes_to_show, s_vals):
        legend_ax.scatter([x_b], [y_cur], s=s_val, color='grey',
                          edgecolor='black', lw=.6, transform=legend_ax.transAxes,
                          clip_on=False)
        legend_ax.text(x_l, y_cur, f"{pct}%", ha='left', va='center',
                       fontsize=label_fs, transform=legend_ax.transAxes)
        y_cur -= row_step  # Space proportional to true bubble size
    # MODIFICATION END

    fig.tight_layout()
    fig.savefig(png, dpi=300, bbox_inches="tight", pad_inches=0.35)
    plt.close(fig);
    print("✓", png)

# ─── GO aspect plotting helpers ─────────────────────────────────
ASPECT_ORDER = [
    "Biological Process",
    "Molecular Function",
    "Cellular Component",
    "Kegg Pathway",
]


def _prepare_top_terms(df, skip_set=""):
    work = df.copy()
    if skip_set:
        work = work[~work["Term"].map(pretty).str.lower().isin(skip_set)]

    ordered_groups = list(GROUPS) + [
        g for g in work["Group"].dropna().unique().tolist() if g not in GROUPS
    ]
    keep = [work[work["Group"] == g].nsmallest(TOP_N, "Pvalue") for g in ordered_groups]
    keep = [k for k in keep if not k.empty]
    if not keep:
        return pd.DataFrame(columns=work.columns)

    top = pd.concat(keep, ignore_index=True)
    top["pretty"] = top["Term"].map(pretty)
    top = (
        top.assign(_group_sort=top["Group"].map(str))
        .sort_values(["_group_sort", "Pvalue"])
        .drop_duplicates("pretty", keep="first")
        .head(MAX_ROWS)
        .drop(columns="_group_sort")
        .copy()
    )
    top["label"] = top["pretty"]
    return top


def _ordered_groups_for_plot(top):
    counts = top["Group"].value_counts()
    min_p = top.groupby("Group")["Pvalue"].min()
    return sorted(
        counts.index.tolist(),
        key=lambda g: (-counts[g], min_p[g], str(g)),
    )


def _style_legend_box(leg):
    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(0.9)
    leg.get_frame().set_facecolor("white")


def _inches_to_axes_frac(fig, ax, inches):
    """Convert physical inches to vertical axes fraction."""
    ax_h_px = ax.get_window_extent().height
    if ax_h_px <= 0: return 0.0
    ax_h_in = ax_h_px / fig.dpi
    return inches / ax_h_in

def _size_legend_layout(fig, legend_ax, sizes_to_show, scale):
    """Compute legend marker areas and absolute spacing in axes-fraction units."""
    # Base layout constants in physical inches
    BASE_ROW_STEP_INCH = 0.25
    MIN_TITLE_GAP_INCH = 0.30
    
    # Calculate marker sizes
    s_vals = np.sqrt(np.asarray(sizes_to_show, dtype=float)) * scale
    
    # Calculate max marker diameter in inches
    # s is area in points^2. diameter in points = sqrt(s). 72 points = 1 inch.
    diams_pts = np.sqrt(s_vals)
    max_diam_pts = diams_pts.max() if len(diams_pts) > 0 else 0
    max_diam_inch = max_diam_pts / 72.0
    
    # Row step: larger of base step or marker diameter + padding
    row_step_inch = max(BASE_ROW_STEP_INCH, max_diam_inch * 1.2)
    
    # Title gap: larger of min gap or half marker + padding
    title_gap_inch = max(MIN_TITLE_GAP_INCH, max_diam_inch * 0.6 + 0.15)
    
    # Convert to axes fraction
    row_step_ax = _inches_to_axes_frac(fig, legend_ax, row_step_inch)
    title_gap_ax = _inches_to_axes_frac(fig, legend_ax, title_gap_inch)
    
    return s_vals, row_step_ax, title_gap_ax


def _draw_go_legends(legend_ax, fig, top, present_groups, scale, fontsize):
    if not present_groups:
        return

    legend_ax.axis("off")
    leg_fs = fontsize * 0.72
    title_fs = leg_fs * 1.20
    x_center = 0.5
    # Keep a small headroom margin so tight bounding never clips top of legend.
    next_y = 0.97

    # Keep only the % genes size legend for GO mode.
    # Use deterministic geometry (instead of artist extents) so circles always stay inside the box.
    sizes_to_show = [100, 75, 50, 25, 10, 5]
    box_w = 0.52
    box_left = x_center - box_w / 2
    box_right = x_center + box_w / 2
    # Layout constants in inches
    TOP_PAD_INCH = 0.2
    BOTTOM_PAD_INCH = 0.15
    
    box_top = next_y
    fig.canvas.draw()
    while True:
        s_vals, row_step, title_gap = _size_legend_layout(
            fig, legend_ax, sizes_to_show, scale
        )
        
        # dynamic padding in axes units
        top_pad = _inches_to_axes_frac(fig, legend_ax, TOP_PAD_INCH)
        bottom_pad = _inches_to_axes_frac(fig, legend_ax, BOTTOM_PAD_INCH)
        
        needed_h = top_pad + title_gap + (max(1, len(sizes_to_show)) - 1) * row_step + bottom_pad
        box_bottom = box_top - needed_h
        if box_bottom >= 0.02 or len(sizes_to_show) <= 2:
            break
        sizes_to_show = sizes_to_show[:-1]

    size_box = FancyBboxPatch(
        (box_left, box_bottom),
        box_w,
        box_top - box_bottom,
        boxstyle="round,pad=0.01",
        transform=legend_ax.transAxes,
        facecolor="white",
        edgecolor="black",
        linewidth=0.9,
        zorder=0,
    )
    legend_ax.add_patch(size_box)

    title_y = box_top - top_pad
    legend_ax.text(
        x_center,
        title_y,
        "% Genes",
        ha="center",
        va="center",
        fontsize=title_fs,
        weight="bold",
        transform=legend_ax.transAxes,
        clip_on=False,
        zorder=2,
    )

    size_box.set_height(box_top - box_bottom)
    size_box.set_y(box_bottom)

    y_cur = title_y - title_gap
    x_b = box_left + 0.13
    x_l = box_left + 0.27
    for pct, s_val in zip(sizes_to_show, s_vals):
        legend_ax.scatter(
            [x_b],
            [y_cur],
            s=float(s_val),
            color="grey",
            edgecolor="black",
            lw=0.6,
            transform=legend_ax.transAxes,
            clip_on=False,
            zorder=2,
        )
        legend_ax.text(
            x_l, y_cur, f"{pct}%",
            ha="left", va="center",
            fontsize=leg_fs * 0.95,
            transform=legend_ax.transAxes,
            zorder=2,
        )
        y_cur -= row_step


def save_aspect_plots_separately(df, out_dir, *, fontsize=16, spacing=10000, skip_set="", plot_width=3):
    """Save one plot per GO aspect using the shared aspect_plot pipeline."""
    out_dir.mkdir(parents=True, exist_ok=True)
    top = _prepare_top_terms(df, skip_set=skip_set)
    if top.empty:
        return

    for aspect in ASPECT_ORDER:
        if not (top["Aspect"] == aspect).any():
            continue
        aspect_plot(
            df,
            out_dir / f"{aspect}.png",
            fontsize=fontsize,
            spacing=spacing,
            skip_set=skip_set,
            plot_width=plot_width,
            aspect_filter=aspect,
        )


def aspect_plot(
    df,
    png,
    *,
    fontsize=16,
    spacing=10000,
    skip_set="",
    plot_width=3,
    aspect_filter=None,
    title=None,
):
    """Bubble plot for all GO aspects or a single aspect when aspect_filter is set."""
    top = _prepare_top_terms(df, skip_set=skip_set)
    if top.empty:
        return

    if aspect_filter is not None:
        top = top[top["Aspect"] == aspect_filter].copy()
        if top.empty:
            return

    top = top.reset_index(drop=True)
    ordered_groups = _ordered_groups_for_plot(top)
    row_step = max(0.8, fontsize / 12 * 0.35)
    group_gap = row_step * 0.8

    y_map = {}
    cur_y = 0.0
    for group_name in reversed(ordered_groups):
        grp_rows = top[top["Group"] == group_name].sort_values("Pvalue", ascending=False)
        for idx in grp_rows.index:
            y_map[idx] = cur_y
            cur_y += row_step
        cur_y += group_gap
    top["y"] = top.index.map(y_map)

    present = [g for g in ordered_groups if g in top["Group"].unique()]
    max_y = top["y"].max() if len(top) else 0
    text_edges = top["-log10P"] + X_PAD + 0.55 * top["label"].str.len()
    max_content_x = text_edges.max()
    for g in present:
        sub = top[top["Group"] == g]
        if sub.empty:
            continue
        group_label = str(GROUP2LABEL.get(g, g))
        bracket_x = (sub["-log10P"] + X_PAD + 0.55 * sub["label"].str.len()).max() + 0.45
        group_label_edge = bracket_x + 0.60 + 0.45 * len(group_label)
        max_content_x = max(max_content_x, group_label_edge)
    xlim = max(5, math.ceil(max_content_x + 0.75))

    n_present = len(present)
    if n_present > 30:
        legend_ratio = 0.75
    elif n_present > 15:
        legend_ratio = 0.6
    else:
        legend_ratio = 0.5

    fig_h = max(6, min(36, len(top) * row_step * 0.40 + 3))
    fig_w = min(24, max(12, xlim / max(1.0, plot_width) + 6))
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, legend_ratio], wspace=0.10)
    ax = fig.add_subplot(gs[0, 0])
    legend_ax = fig.add_subplot(gs[0, 1])

    scale = 650 / np.sqrt(top["%"].max() or 1)

    for g in ordered_groups:
        sub = top[top["Group"] == g]
        if sub.empty:
            continue
        face = GROUP2HEX.get(g, "#7f7f7f")
        if is_no_broad_group(g):
            face = "white"
        ax.scatter(
            sub["-log10P"], sub["y"],
            s=np.sqrt(sub["%"]) * scale, color=face,
            edgecolor="black", lw=0.6, zorder=3
        )

    for _, row in top.iterrows():
        ax.text(
            row["-log10P"] + X_PAD, row["y"], row["label"],
            ha="left", va="center", fontsize=fontsize * 0.75,
            color="#222222",
            zorder=4, path_effects=TEXT_OUTLINE
        )

    # Broad-category brackets and labels (muted so GO terms remain primary)
    bracket_alpha = 0.45
    bracket_text_alpha = 0.55
    for g in present:
        if is_no_broad_group(g):
            continue
        sub = top[top["Group"] == g]
        if sub.empty:
            continue
        base = GROUP2HEX.get(g, "#7f7f7f")
        y0, y1 = sub["y"].min() - 0.35, sub["y"].max() + 0.35
        mid_y = (y0 + y1) / 2
        label_edge = (sub["-log10P"] + X_PAD + 0.55 * sub["label"].str.len()).max()
        bracket_x = label_edge + 0.45
        ax.plot([bracket_x, bracket_x], [y0, y1], color=base, lw=1.15, alpha=bracket_alpha, zorder=2)
        ax.plot([bracket_x - 0.75, bracket_x], [y0, y0], color=base, lw=1.15, alpha=bracket_alpha, zorder=2)
        ax.plot([bracket_x - 0.75, bracket_x], [y1, y1], color=base, lw=1.15, alpha=bracket_alpha, zorder=2)
        ax.text(
            bracket_x + 0.60,
            mid_y,
            str(GROUP2LABEL.get(g, g)),
            ha="left",
            va="center",
            fontsize=fontsize * 0.70,
            weight="bold",
            color=base,
            alpha=bracket_text_alpha,
            zorder=2,
        )

    ax.set_xlim(0, xlim)
    ax.set_ylim(-0.5, max_y + 0.5)
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    ax.tick_params(axis="x", labelsize=fontsize * 0.75, direction="out", length=4, width=0.8)
    ax.grid(axis="x", ls="--", alpha=0.12)
    ax.set_xlabel("-log10(P)", fontsize=fontsize * 0.9)
    if title:
        ax.set_title(title, weight="bold", fontsize=fontsize * 0.95)

    _draw_go_legends(legend_ax, fig, top, present, scale, fontsize)

    fig.savefig(png, dpi=300, bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)
    print("✓", png)

# ─── plotting (vertical) ────────────────────────────────────────
def vertical_plot(w_df, s_df, png, *, fontsize=16, skip_set="", plot_width=3):
    """
    Creates a vertically stacked bubble plot for wake vs. sleep data with aligned brackets.
    """

    # --- Helper to get the data that will be plotted ---
    def _get_top_df(df):
        if skip_set:
            df = df[~df["Term"].map(pretty).str.lower().isin(skip_set)]
        keep = [df[df["Group"] == g].nsmallest(TOP_N, "Pvalue") for g in GROUPS]
        top = pd.concat(keep, ignore_index=True)
        if top.empty: return pd.DataFrame(columns=df.columns)
        top["pretty"] = top["Term"].map(pretty)
        top = (top.sort_values(["Group", "Pvalue"])
               .drop_duplicates("pretty", keep="first")
               .head(MAX_ROWS))
        top["label"] = top["pretty"]
        return top

    w_top = _get_top_df(w_df)
    s_top = _get_top_df(s_df)
    all_top = pd.concat([w_top, s_top])
    if all_top.empty: return

    # --- Calculate a shared X-axis limit and a fixed bracket position ---
    # 1. Find the maximum extent of the bubble labels
    max_text_edge = (all_top["-log10P"] + X_PAD + LABEL_FACTOR * all_top["label"].str.len()).max()

    # 2. Define a fixed amount of space for the brackets and their labels
    bracket_zone_width =20
    final_xlim = math.ceil((max_text_edge + bracket_zone_width) / 5.0) * 5 # Round up to nearest 5

    # 3. Set a single, fixed x-position for all brackets
    bracket_xb = final_xlim - (bracket_zone_width * 0.8)


    # --- Panel-drawing helper function ---
    def _draw_panel(ax, top_df, title):
        if top_df.empty: return

        y_map = pack_y_grouped(top_df, max(.8, fontsize / 12 * .35))
        if not y_map: return
        top_df["y"] = top_df["pretty"].map(y_map)
        max_y = max(y_map.values())
        scale = 650 / np.sqrt(top_df["%"].max() or 1)

        for g in GROUPS:
            sub = top_df[top_df["Group"] == g]
            if sub.empty: continue
            base = GROUP2HEX[g]
            
            ax.scatter(sub["-log10P"], sub["y"], s=np.sqrt(sub["%"]) * scale, color=base,
                       edgecolor='black', lw=.6, zorder=3)

            for x, y, lbl, aspect in zip(sub["-log10P"], sub["y"], sub["label"], sub.get("Aspect", pd.Series([""] * len(sub)))):
                ax.text(x + X_PAD, y, lbl, ha='left', va='center', fontsize=fontsize * .8,
                        color=aspect_text_color(aspect), zorder=4, path_effects=TEXT_OUTLINE)

            # Use the single, pre-calculated, shared bracket position for all groups
            y0, y1 = sub["y"].min() - 0.4, sub["y"].max() + 0.4
            mid = (y0 + y1) / 2
            ax.plot([bracket_xb, bracket_xb], [y0, y1], color=base, lw=1.2)
            ax.plot([bracket_xb - 1, bracket_xb], [y0, y0], color=base, lw=1.2)
            ax.plot([bracket_xb - 1, bracket_xb], [y1, y1], color=base, lw=1.2)
            ax.text(bracket_xb + 0.8, mid, GROUP2LABEL[g], ha='left', va='center',
                    fontsize=fontsize * .85, color=base, weight='bold')

        ax.set_ylim(-1, max_y + 1)
        ax.set_yticks([])
        ax.set_title(title, weight='bold', fontsize=fontsize)

    # --- Main plotting logic ---
    fig_w, fig_h = 24, 18
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 0.25], wspace=0.15, hspace=0.2)
    ax_w = fig.add_subplot(gs[0, 0])
    ax_s = fig.add_subplot(gs[1, 0], sharex=ax_w)
    legend_ax = fig.add_subplot(gs[:, 1])

    _draw_panel(ax_w, w_top, "Wake Genes GO Terms")
    _draw_panel(ax_s, s_top, "Sleep Genes GO Terms")

    # Configure shared X-axis
    ax_s.set_xlim(0, final_xlim)  # sharex means this sets both
    ax_s.set_xlabel("Negative log p-value", fontsize=fontsize + 6, weight='bold')
    ax_s.tick_params(axis='x', labelsize= 24)
    plt.setp(ax_w.get_xticklabels(), visible=False)
    ax_w.grid(axis='x', ls='--', alpha=.15)
    ax_s.grid(axis='x', ls='--', alpha=.15)

    # --- Legend (unchanged) ---
    legend_ax.axis("off")
    present_groups = all_top["Group"].unique()
    present = [g for g in GROUPS if g in present_groups][::-1]
    if not present:
        fig.tight_layout()
        fig.savefig(png, dpi=300, bbox_inches="tight", pad_inches=0.35);
        plt.close(fig);
        return

    handles = [plt.Line2D([0], [0], marker='o', ls='', color=GROUP2HEX[g], markersize=8) for g in present]
    leg = legend_ax.legend(handles, [GROUP2LABEL[g] for g in present], title="Group", loc="upper left",
                           frameon=False, fontsize=fontsize * .75, title_fontsize=fontsize * .8)
    legend_ax.add_artist(leg)
    fig.canvas.draw()
    y_cur = leg.get_window_extent().transformed(legend_ax.transAxes.inverted()).y0
    t_fs, label_fs = fontsize * .8, fontsize * .75

    y_cur -= 0.10
    legend_ax.text(.5, y_cur, "% Genes", ha='center', va='bottom', fontsize=t_fs, transform=legend_ax.transAxes)
    y_cur -= 0.04
    sizes_to_show = [100, 75, 50, 25, 10, 5]
    x_b, x_l = .35, .65

    # Use combined data for scale calculation
    max_pct = all_top['%'].max()
    scale = 650 / np.sqrt(max_pct or 1)
    fig.canvas.draw()
    bottom_limit = 0.02
    while True:
        s_vals, row_step, title_gap = _size_legend_layout(fig, legend_ax, sizes_to_show, scale)
        first_y = y_cur - title_gap
        last_y = first_y - (len(sizes_to_show) - 1) * row_step
        if last_y >= bottom_limit or len(sizes_to_show) <= 2:
            break
        sizes_to_show = sizes_to_show[:-1]
    y_cur = first_y
    for pct, s_val in zip(sizes_to_show, s_vals):
        legend_ax.scatter([x_b], [y_cur], s=s_val, color='grey', edgecolor='black', lw=.6,
                          transform=legend_ax.transAxes, clip_on=False)
        legend_ax.text(x_l, y_cur, f"{pct}%", ha='left', va='center', fontsize=label_fs, transform=legend_ax.transAxes)
        y_cur -= row_step

    fig.tight_layout()
    fig.savefig(png, dpi=300, bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)
    print("✓", png)

# ─── CLI driver ─────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser()
    p.add_argument("indir")
    p.add_argument("--fontsize", type=int, default=16)
    p.add_argument("--spacing",  type=int, default=10000)
    p.add_argument("--vertical", action="store_true",
                   help="stack sleep & wake panels vertically instead of side-by-side")
    p.add_argument("--width", type = float, default = 3)
    p.add_argument("--keep-uncat", action="store_true",
                   help="also plot GO terms that do not match any category")
    p.add_argument("--by_GO", action="store_true",
                   help="ignore keyword categories and colour bubbles by GO aspect (BP/MF/CC) or KEGG pathway")
    p.add_argument("--go-map-file", type=pathlib.Path, default=CLUSTER_FILE,
                   help="path to GO-term mapping file used with --by_GO")
    p.add_argument("--split-terms", action="store_true",
                   help="generate one PNG per category (including 'either' if present) in a folder named after the file")
    a = p.parse_args()

    indir=pathlib.Path(a.indir)
    out=indir/"Plots"; out.mkdir(exist_ok=True)

    def _norm_group(g):
        # keep numeric groups as python ints (not numpy scalar wrappers)
        if isinstance(g, np.integer):
            return int(g)
        if isinstance(g, np.floating) and float(g).is_integer():
            return int(g)
        return g

    def _ensure_group_colored(_g):
        return

    # Initialize ClusterAssigner
    cluster_assigner = None
    if a.by_GO:
        cluster_assigner = ClusterAssigner(a.go_map_file, unknown_policy="unique")
        GROUPS.clear()
        GROUPS.extend(cluster_assigner.known_groups)
        if "Unassigned" not in GROUPS: GROUPS.append("Unassigned")
        
        GROUP2KW.clear()
        GROUP2HEX.clear()
        GROUP2LABEL.clear()
        
        # Generate colors
        from matplotlib.colors import to_hex
        
        n_groups = len(GROUPS)
        colors = []
        
        # Strategy for colors:
        # <= 20: tab20
        # <= 40: tab20 + tab20b
        # <= 60: tab20 + tab20b + tab20c
        # > 60: nipy_spectral (high contrast rainbow)
        
        if n_groups <= 20:
            cmap = plt.get_cmap("tab20")
            colors = [to_hex(cmap(i)) for i in range(20)]
        elif n_groups <= 40:
            c1 = [to_hex(plt.get_cmap("tab20")(i)) for i in range(20)]
            c2 = [to_hex(plt.get_cmap("tab20b")(i)) for i in range(20)]
            colors = c1 + c2
        elif n_groups <= 60:
            c1 = [to_hex(plt.get_cmap("tab20")(i)) for i in range(20)]
            c2 = [to_hex(plt.get_cmap("tab20b")(i)) for i in range(20)]
            c3 = [to_hex(plt.get_cmap("tab20c")(i)) for i in range(20)]
            colors = c1 + c2 + c3
        else:
            cmap = plt.get_cmap("nipy_spectral")
            # Sample evenly
            colors = [to_hex(cmap(i/n_groups)) for i in range(n_groups)]
            
        for i, g in enumerate(GROUPS):
            GROUP2HEX[g] = colors[i % len(colors)]
            GROUP2LABEL[g] = g

        def _ensure_group_colored(g):
            # Add to GROUPS if missing
            g = _norm_group(g)
            if is_no_broad_group(g):
                return
            if g not in GROUPS:
                GROUPS.append(g)
            # Assign a color/label if missing
            if g not in GROUP2HEX:
                idx = len(GROUP2HEX)  # next color index
                GROUP2HEX[g] = colors[idx % len(colors)]
            if g not in GROUP2LABEL:
                GROUP2LABEL[g] = g


    # optional skip list
    skip_file=indir/"skip_terms.xlsx"
    skip_set=set()
    # if skip_file.is_file():
    #     skip_terms=pd.read_excel(skip_file)["Keywords"].str.strip().str.lower()
    #     skip_set=set(skip_terms)

    xlsxs=list(indir.glob("*.xlsx"))
    sleep=next((f for f in xlsxs if re.search("sleep",f.stem,re.I)),None)
    wake =next((f for f in xlsxs if re.search("wake", f.stem,re.I)),None)

    for f in xlsxs:
        # if f in (sleep,wake): continue
        active_cluster_assigner = cluster_assigner
        if a.by_GO:
            active_cluster_assigner = ClusterAssigner(a.go_map_file, unknown_policy="unique", source_name=f.stem)
            validate_all_terms_mapped([f], active_cluster_assigner)
        df = read_csv(f, keep_uncat=a.keep_uncat, aspect_mode=a.by_GO, cluster_assigner=active_cluster_assigner)
        if a.by_GO:
            for g in df["Group"].dropna().unique():
                _ensure_group_colored(g)
        if a.split_terms:
            # ... (existing code for split terms) ...
            # Create folder named after file in GO_Analysis/Plots folder within input directory
            split_dir = indir / "Plots" / f.stem
            split_dir.mkdir(parents=True, exist_ok=True)
            
            # Get all groups present in the data (including "either" if it exists as a category)
            present_groups = sorted(df["Group"].dropna().unique().tolist())
            
            # Generate one plot per category/group
            for group in present_groups:
                # Filter dataframe to only this group
                group_df = df[df["Group"] == group].copy()
                if group_df.empty:
                    continue
                
                # Clean group name for filename (replace spaces/special chars)
                safe_group_name = re.sub(r'[^\w\s-]', '', str(group)).strip()
                
                # Apply TOP_N filtering first (same as bubble_plot does internally)
                # This ensures we only split the terms that would actually be plotted
                filtered_group_df = group_df.nsmallest(min(TOP_N, len(group_df)), "Pvalue").copy()
                filtered_group_df["pretty"] = filtered_group_df["Term"].map(pretty)
                filtered_group_df = (filtered_group_df.sort_values("Pvalue")
                                     .drop_duplicates("pretty", keep="first"))
                
                # Check if we need to split across multiple plots
                num_terms = len(filtered_group_df)
                if num_terms <= MAX_TERMS_PER_PLOT:
                    # Single plot
                    output_file = split_dir / f"{safe_group_name}.png"
                    bubble_plot(filtered_group_df, output_file,
                               fontsize=a.fontsize, spacing=a.spacing, skip_set=skip_set, plot_width=a.width)
                else:
                    # Split into multiple plots
                    num_plots = math.ceil(num_terms / MAX_TERMS_PER_PLOT)
                    for plot_num in range(1, num_plots + 1):
                        start_idx = (plot_num - 1) * MAX_TERMS_PER_PLOT
                        end_idx = plot_num * MAX_TERMS_PER_PLOT
                        chunk_df = filtered_group_df.iloc[start_idx:end_idx].copy()
                        
                        if chunk_df.empty:
                            continue
                        
                        output_file = split_dir / f"{safe_group_name}{plot_num}.png"
                        bubble_plot(chunk_df, output_file,
                                   fontsize=a.fontsize, spacing=a.spacing, skip_set=skip_set, plot_width=a.width)
        else:
            if a.by_GO:
                # Use the new separate files saving function
                # Output directory: Plots/FilenameStem/
                # Note: 'out' is currently Plots/
                file_stem_dir = out / f.stem
                
                # 1. Save separate aspect plots
                save_aspect_plots_separately(df, file_stem_dir,
                           fontsize=a.fontsize, spacing=a.spacing, skip_set=skip_set, plot_width=a.width)
                
                # 2. Save combined plot (original aspect_plot)
                # Save it in the subdirectory as "All_Aspects.png"
                combined_png_path = file_stem_dir / "All_Aspects.png"
                aspect_plot(df, combined_png_path,
                           fontsize=a.fontsize, spacing=a.spacing, skip_set=skip_set, plot_width=a.width)
            else:
                # Original behavior: single PNG per file
                bubble_plot(df, out/f"{f.stem}.png",
                           fontsize=a.fontsize, spacing=a.spacing, skip_set=skip_set, plot_width=a.width)

if __name__=="__main__":
    main()
