#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
from pathlib import Path


INPUT_FILE = "data/processed/precision_analysis.csv"
OUTPUT_DIR = "results/figures"
OUTPUT_BASE = "Figure3"

# Coverage range colors (green = high coverage, red = zero)
COVERAGE_COLORS = {
    "90-100": "#00441B",
    "80-90": "#1B7837",
    "50-80": "#A6DBA0",
    "10-50": "#D9F0D3",
    "<10": "#FEE08B",
    "0": "#D73027"
}

COVERAGE_RANK = {"0": 0, "<10": 1, "10-50": 2, "50-80": 3, "80-90": 4, "90-100": 5}

# Run group labels: Combined Intra-pool and Inter-pool into "Intra-assay"
RUN_GROUP_LABELS = {
    "Within run, within pool": "Intra-assay",
    "Within run, between pool": "Intra-assay",
    "Between runs": "Inter-assay"
}

# Plot styling
PLOT_CONFIG = {
    "font_family": "Arial",
    "font_size": 10,
    "title_size": 18,
    "label_size": 18,
    "tick_size": 14,
    "cell_size_inches": 0.6,
    "dpi": 300,
    "vertical_gap": 0.3,      # Gap between subplots in inches
    "top_margin": 0.8,        # Top margin in inches (for Ct labels)
    "bottom_margin": 0.8,     # Bottom margin in inches (for x-axis label)
    "left_margin": 0.8,       # Left margin in inches (for y-axis label)
    "right_margin": 0.5,      # Right margin in inches
    "run_group_title_height": 0.3  # Height reserved for run group titles
}

def load_and_process_data(csv_file):

    df = pd.read_csv(csv_file)
    
    # Convert coverage to categorical ranges
    def get_coverage_range(coverage):
        if coverage == 0:
            return "0"
        elif coverage < 10:
            return "<10"
        elif coverage < 50:
            return "10-50"
        elif coverage < 80:
            return "50-80"
        elif coverage < 90:
            return "80-90"
        else:
            return "90-100"
    
    # Categorize by Ct value
    def get_ct_category(ct_value):
        return "Ct <30" if ct_value < 30 else "Ct >30"
    
    # Classify run groups based on Run/Pool combinations per sample
    def classify_run_groups(df):
        sample_classifications = {}
        for sample_id in df['Sample-id'].unique():
            sample_data = df[df['Sample-id'] == sample_id]
            combos = sample_data[['Run', 'Pool']].drop_duplicates()
            
            if len(combos) > 1:
                if len(combos['Run'].unique()) > 1:
                    sample_classifications[sample_id] = "Between runs"
                else:
                    sample_classifications[sample_id] = "Within run, between pool"
            else:
                sample_classifications[sample_id] = "Within run, within pool"
        
        return sample_classifications
    
    # Apply transformations
    df['Coverage_range'] = df['Coverage'].apply(get_coverage_range)
    df['CT_category'] = df['CT-value'].apply(get_ct_category)
    df['Run_group'] = df['Sample-id'].map(classify_run_groups(df))
    df['Run_group_relabeled'] = df['Run_group'].map(RUN_GROUP_LABELS)
    df['coverage_rank'] = df['Coverage_range'].map(COVERAGE_RANK).fillna(0).astype(int)
    df['Virus'] = df['Virus'].replace({"HRV": "Rhinovirus"})
    
    # Keep highest coverage for duplicate entries
    idx = df.groupby(['Virus', 'Sample-id', 'CT_category', 'Workflow', 'Run_group'])['coverage_rank'].idxmax()
    df = df.loc[idx].reset_index(drop=True)
    
    return df


def create_triangle_vertices(col, row, position="top_right"):

    if position == "top_right":
        return [(col + 1, row + 1), (col, row + 1), (col + 1, row)]
    else:
        return [(col, row), (col + 1, row), (col, row + 1)]


def create_precision_plot(df, output_dir, output_base):

    plt.rcParams.update({
        "font.family": PLOT_CONFIG["font_family"],
        "font.size": PLOT_CONFIG["font_size"],
        "axes.titlesize": PLOT_CONFIG["title_size"],
        "axes.labelsize": PLOT_CONFIG["label_size"],
        "xtick.labelsize": PLOT_CONFIG["tick_size"],
        "ytick.labelsize": PLOT_CONFIG["tick_size"],
        "axes.linewidth": 1.5,
        "axes.titleweight": "bold",
        "axes.labelweight": "bold",
    })
    
    # Define run groups: Only Intra-assay and Inter-assay
    run_groups = ["Intra-assay", "Inter-assay"]
    
    # Get all viruses
    viruses = sorted(df['Virus'].unique())
    n_rows = len(viruses)
    
    # Determine samples per panel
    panel_samples = {}
    max_samples = 0
    for run_group in run_groups:
        subset = df[df['Run_group_relabeled'] == run_group]
        samples = sorted(subset['Sample-id'].unique())
        panel_samples[run_group] = samples
        max_samples = max(max_samples, len(samples))
    
    # Calculate figure dimensions for perfect squares
    cell_width_inches = PLOT_CONFIG["cell_size_inches"]
    cell_height_inches = cell_width_inches  # Square cells
    
    # Calculate width for each panel
    panel_widths = [len(panel_samples[g]) * cell_width_inches for g in run_groups]
    total_width = sum(panel_widths) + PLOT_CONFIG["left_margin"] + PLOT_CONFIG["right_margin"]
    
    # Calculate height for each subplot
    subplot_height = n_rows * cell_height_inches + PLOT_CONFIG["run_group_title_height"]
    total_height = (2 * subplot_height + 
                   PLOT_CONFIG["vertical_gap"] + 
                   PLOT_CONFIG["top_margin"] + 
                   PLOT_CONFIG["bottom_margin"])
    
    # Create figure
    fig = plt.figure(figsize=(total_width, total_height))
    
    # Calculate normalized width ratios for GridSpec
    width_ratios = [len(panel_samples.get(g, [])) for g in run_groups]
    width_ratios = [w if w > 0 else 1 for w in width_ratios]
    
    # Calculate normalized height ratios (subplot heights)
    subplot_height_ratio = subplot_height
    height_ratios = [subplot_height_ratio, subplot_height_ratio]
    
    # Create GridSpec
    gs = plt.GridSpec(2, len(run_groups), figure=fig,
                     hspace=PLOT_CONFIG["vertical_gap"] / subplot_height,
                     wspace=0.05,
                     height_ratios=height_ratios,
                     width_ratios=width_ratios,
                     left=PLOT_CONFIG["left_margin"] / total_width,
                     right=1 - (PLOT_CONFIG["right_margin"] / total_width),
                     top=1 - (PLOT_CONFIG["top_margin"] / total_height),
                     bottom=PLOT_CONFIG["bottom_margin"] / total_height)
    
    # Create axes
    axes_top = [fig.add_subplot(gs[0, i]) for i in range(len(run_groups))]
    axes_bottom = [fig.add_subplot(gs[1, i]) for i in range(len(run_groups))]
    
    # Function to plot a single panel
    def plot_panel(ax, df_filtered, sample_ids, viruses, n_rows, n_samples, show_x_labels=False):
        if not sample_ids or len(sample_ids) == 0:
            ax.set_visible(False)
            return
        
        # Light grey grid borders (behind everything - zorder=1)
        for col_idx in range(n_samples):
            rect = patches.Rectangle((col_idx, 0), 1, n_rows, fill=False,
                                    edgecolor='#E6E6E6', linewidth=2.0, zorder=1)
            ax.add_patch(rect)
        
        for row_idx in range(n_rows):
            rect = patches.Rectangle(
                (0, row_idx), n_samples, 1,
                fill=False,
                edgecolor='#E6E6E6',
                linewidth=2.0,
                zorder=1
            )
            ax.add_patch(rect)
        
        # Draw cells
        for row_idx, virus in enumerate(viruses):
            for col_idx, sample_id in enumerate(sample_ids):
                cell = df_filtered[(df_filtered['Virus'] == virus) & 
                                  (df_filtered['Sample-id'] == sample_id)]
                
                has_data = False
                
                # Check if both workflows have coverage 90-100
                r_id = cell[cell['Workflow'] == 'ID']
                r_idd = cell[cell['Workflow'] == 'ID-D']
                both_darkest_green = (not r_id.empty and r_id.iloc[0]['Coverage_range'] == "90-100" and 
                                      not r_idd.empty and r_idd.iloc[0]['Coverage_range'] == "90-100")
                
                # Set diagonal line color based on condition
                diagonal_color = 'white' if both_darkest_green else 'black'
                
                # ID workflow (top-right triangle)
                if not r_id.empty:
                    color = COVERAGE_COLORS.get(r_id.iloc[0]['Coverage_range'], "#FFFFFF")
                    poly = patches.Polygon(
                        create_triangle_vertices(col_idx, row_idx, "top_right"),
                        closed=True, facecolor=color, edgecolor=diagonal_color, linewidth=2.0, zorder=2
                    )
                    ax.add_patch(poly)
                    has_data = True
                
                # ID-D workflow (bottom-left triangle)
                if not r_idd.empty:
                    color = COVERAGE_COLORS.get(r_idd.iloc[0]['Coverage_range'], "#FFFFFF")
                    poly = patches.Polygon(
                        create_triangle_vertices(col_idx, row_idx, "bottom_left"),
                        closed=True, facecolor=color, edgecolor=diagonal_color, linewidth=2.0, zorder=2
                    )
                    ax.add_patch(poly)
                    has_data = True
                
                # Black cell border (on top of grey grid - zorder=3)
                if has_data:
                    rect = patches.Rectangle((col_idx, row_idx), 1, 1, fill=False,
                                            edgecolor='black', linewidth=1.8, zorder=3)
                    ax.add_patch(rect)
        
        # Panel border
        ax.add_patch(patches.Rectangle((0, 0), n_samples, n_rows, fill=False, edgecolor='grey', linewidth=3.0, zorder=4))
        
        # Axes configuration
        ax.set_xlim(0, n_samples)
        ax.set_ylim(0, n_rows)
        ax.set_xlim(0, n_samples)
        ax.set_ylim(0, n_rows)
        
        # Set x-axis ticks (only for bottom subplot)
        if show_x_labels:
            ax.set_xticks([i + 0.5 for i in range(n_samples)])
            ax.set_xticklabels(sample_ids, ha='right', fontsize=12, rotation=45, rotation_mode='anchor')
            ax.set_xlabel('')
        else:
            ax.set_xticks([i + 0.5 for i in range(n_samples)])
            ax.set_xticklabels([])
        
        # Set y-axis ticks
        ax.set_yticks([j + 0.5 for j in range(n_rows)])
        ax.set_yticklabels(viruses, fontsize=14)
        
        ax.invert_yaxis()
        ax.tick_params(axis='both', length=0)
    
    # Plot top panels (Ct <30)
    for idx_ax, run_group in enumerate(run_groups):
        ax = axes_top[idx_ax]
        subset = df[(df['Run_group_relabeled'] == run_group) & (df['CT_category'] == "Ct <30")]
        sample_ids = panel_samples.get(run_group, [])
        n_samples = len(sample_ids)
        
        plot_panel(ax, subset, sample_ids, viruses, n_rows, n_samples, show_x_labels=False)
        
        # Add run group title at the top of each panel (only for top subplot)
        if n_samples > 0:
            ax.set_title(run_group, fontsize=14, pad=10, fontweight='bold')
        
        # Remove y-axis labels for non-first columns
        if idx_ax != 0:
            ax.set_yticklabels([])
            ax.set_ylabel('')
    
    # Plot bottom panels (Ct >30)
    for idx_ax, run_group in enumerate(run_groups):
        ax = axes_bottom[idx_ax]
        subset = df[(df['Run_group_relabeled'] == run_group) & (df['CT_category'] == "Ct >30")]
        sample_ids = panel_samples.get(run_group, [])
        n_samples = len(sample_ids)
        
        # Show x-axis labels only for bottom subplot
        plot_panel(ax, subset, sample_ids, viruses, n_rows, n_samples, show_x_labels=True)
        
        # Remove y-axis labels for non-first columns
        if idx_ax != 0:
            ax.set_yticklabels([])
            ax.set_ylabel('')
    
    # Add single shared x-axis label at the bottom
    fig.text(0.5, -0.05, 'Sample ID', ha='center', fontsize=18, fontfamily='Arial', fontweight='bold')
    
    # Add single shared y-axis label on the left
    fig.text(-0.05, 0.5, 'Virus', ha='center', va='center', rotation=90, fontsize=18, fontfamily='Arial', fontweight='bold')
    
    # Add legend
    legend_elements = [
        Line2D([0], [0], marker='s', color='w', label=k,
               markerfacecolor=v, markersize=16, markeredgecolor='grey', 
               linewidth=1.5, markeredgewidth=1.5)
        for k, v in COVERAGE_COLORS.items()
    ]
    
    legend = fig.legend(
        handles=legend_elements,
        title="Genome Coverage Range (%)",
        loc='lower right',
        ncol=len(COVERAGE_COLORS),
        bbox_to_anchor=(0.94, -0.15),
        frameon=False,
        fontsize=20,
        handletextpad=0.8,
        columnspacing=0.3,
        ncols=3
    )
    legend.get_title().set_fontsize(20)
    
    # Save outputs
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    plt.savefig(output_path / f"{output_base}.png", dpi=PLOT_CONFIG["dpi"], bbox_inches="tight")
    plt.savefig(output_path / f"{output_base}.pdf", bbox_inches="tight")
    plt.close()


def create_legend_only(output_dir, filename="figure3_legend.png"):

    fig_legend = plt.figure(figsize=(8, 1))
    
    legend_elements = [
        Line2D([0], [0], marker='s', color='w', label=k,
               markerfacecolor=v, markersize=10, markeredgecolor='grey', linewidth=0)
        for k, v in COVERAGE_COLORS.items()
    ]
    
    legend = fig_legend.legend(
        handles=legend_elements,
        title="Coverage Range (%)",
        loc="center",
        ncol=len(COVERAGE_COLORS),
        frameon=True,
        fancybox=False,
        edgecolor='black',
        fontsize=10
    )
    legend.get_title().set_fontsize(11)
    legend.get_title().set_fontweight('bold')
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    fig_legend.savefig(output_path / filename, dpi=300, bbox_inches='tight',
                       bbox_extra_artists=[legend], pad_inches=0.1)
    plt.close(fig_legend)

def main():
    
    # Validate input file
    if not Path(INPUT_FILE).exists():
        raise FileNotFoundError(
            f"Input file not found: {INPUT_FILE}\n"
            f"Expected columns: Sample-id, Virus, Coverage, CT-value, Workflow, Run, Pool"
        )
    
    df = load_and_process_data(INPUT_FILE)

    create_precision_plot(df, OUTPUT_DIR, OUTPUT_BASE)
    create_legend_only(OUTPUT_DIR)

if __name__ == "__main__":
    main()