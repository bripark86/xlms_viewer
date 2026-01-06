import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import os
import re
import subprocess
from io import StringIO

# Page configuration
st.set_page_config(page_title="2D Data Explorer", layout="wide")

# ===================================================================
# SECTION 1: HELPER FUNCTIONS & CONSTANTS
# ===================================================================

# Amino Acid Coloring Logic
AA_GROUPS = {
    'hydrophobic': ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'P', 'G'],
    'polar': ['S', 'T', 'N', 'Q', 'C'],
    'acidic': ['D', 'E'],
    'basic': ['R', 'K', 'H']
}

def get_aa_class(aa):
    """Get CSS class for amino acid based on its chemical properties."""
    aa_upper = aa.upper()
    for group_name, aa_list in AA_GROUPS.items():
        if aa_upper in aa_list:
            return f"aa-{group_name}"
    return "aa-unknown"

def parse_fasta(lines):
    """Parse FASTA file content and return dictionary of sequences."""
    sequences = {}
    current_header = None
    current_seq = []
    
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if current_header and current_seq:
                sequences[current_header] = "".join(current_seq).upper()
            current_header = line
            current_seq = []
        else:
            if current_header:
                # Remove non-alphabetic characters
                clean_line = re.sub(r'[^a-zA-Z]', '', line)
                if clean_line:
                    current_seq.append(clean_line)
    
    # Handle last sequence
    if current_header and current_seq:
        sequences[current_header] = "".join(current_seq).upper()
    
    return sequences

# Fixed color palette (matching R version)
FIXED_PALETTE = {
    "SMARCA4": "#558B2F", "SMCA4": "#558B2F", "BRG1": "#558B2F",
    "SMARCA2": "#6B8E23", "SMCA2": "#6B8E23", "BRM": "#6B8E23",
    "ARID1A": "#A44A4A", "ARI1A": "#A44A4A",
    "ARID1B": "#CD5C5C", "ARI1B": "#CD5C5C",
    "ARID2": "#F08080",
    "BCL7A": "#483D8B", "BCL7B": "#6A5ACD", "BCL7C": "#7B68EE",
    "SMARCB1": "#4682B4", "SNF5": "#4682B4",
    "SMARCC1": "#5F9EA0", "SMRC1": "#5F9EA0", "BAF155": "#5F9EA0",
    "SMARCC2": "#87CEEB", "SMRC2": "#87CEEB", "BAF170": "#87CEEB",
    "SMARCD1": "#B8860B", "SMRD1": "#B8860B",
    "SMARCD2": "#DAA520", "SMRD2": "#DAA520",
    "SMARCD3": "#FFD700", "SMRD3": "#FFD700",
    "SMARCE1": "#8B4513", "SMCE1": "#8B4513", "BAF57": "#8B4513",
    "SS18": "#A0522D", "SS18L1": "#D2691E",
    "ACTB": "#2E8B57", "ACTL6A": "#3CB371", "ACL6A": "#3CB371", "ACTL6B": "#90EE90",
    "DPF1": "#C71585", "DPF2": "#DB7093", "DPF3": "#FF69B4",
    "PBRM1": "#663399", "PHF10": "#8A2BE2", "BRD7": "#9400D3",
    "BICRA": "#FF4500", "BICRAL": "#FF6347", "BRD9": "#FF7F50"
}

# File info mapping (exact port from R)
FILE_INFO_LIST = {
    "both_cbaf_pbaf": {"stem": "Both_CBAF_and_PBAF_xi_net_between_subunit", "type": "standard"},
    "canonical_between": {"stem": "Canonical_xi_net_between_subunit", "type": "standard"},
    "canonical_intra": {"stem": "Canonical_xi_net_inter_and_intra", "type": "standard"},
    "nucleosome": {"stem": "Nucleosome_and_cBAF_xi_net_between_subunit", "type": "standard"},
    "pbaf": {"stem": "PBAF_xi_net_between_subunit", "type": "standard"},
    "ncbaf_inter1": {"stem": "ncBAF_inter1_links", "type": "standard_shortname"},
    "ncbaf_intra1": {"stem": "ncBAF_intra1_links", "type": "standard_shortname"},
    "cbaf_gr_proxl": {"stem": "cBAF_GR_xlinks-proteins-search-2454-2025-10-09", "type": "standard_shortname"},
    "cbaf_gr_ncp_proxl": {"stem": "cBAF_GR_NCP_xlinks-proteins-search-2515-2025-10-09", "type": "standard_shortname"},
    "ha_dpf2_proxl": {"stem": "HA_DPF2_xlinks-proteins-search-254-2025-10-09", "type": "standard_shortname"},
    "ha_dpf2_ncp_proxl": {"stem": "HA_DPF2_NCP_xlinks-proteins-search-724-2025-10-09", "type": "standard_shortname"},
    "ss18_flag_proxl": {"stem": "SS18_Flag_xlinks-proteins-search-257-2025-10-09", "type": "standard_shortname"},
    "phf10_ha_proxl": {"stem": "PHF10_HA_xlinks-proteins-search-261-2025-10-09", "type": "standard_shortname"},
    "ha_brd7_proxl": {"stem": "HA_BRD7_xlinks-proteins-search-259-2025-10-09", "type": "standard_shortname"}
}

# CSS for sequence viewer
SEQUENCE_VIEWER_CSS = """
<style>
.aa-residue {
    font-family: 'Courier New', Courier, monospace;
    font-weight: bold;
    font-size: 0.8em;
    padding: 1px;
    margin: 1px;
    border: 2px solid transparent;
    border-radius: 1px;
    display: inline-block;
    min-width: 5px;
    text-align: center;
    cursor: default;
    transition: transform 0.1s ease, box-shadow 0.1s ease;
    line-height: 1;
    color: transparent;
    position: relative;
}
.aa-residue:hover {
    transform: scale(1.5);
    color: inherit;
    box-shadow: 0 0 8px rgba(0,0,0,0.5);
    position: relative;
    z-index: 10;
}
.aa-highlight {
    border: 2px solid #fbbf24;
    box-shadow: 0 0 5px #fbbf24;
    transform: scale(1.1);
}
.aa-highlight:hover {
    transform: scale(1.5);
    z-index: 20;
}
.custom-tooltip {
    visibility: hidden;
    opacity: 0;
    background-color: #333;
    color: #fff;
    text-align: center;
    border-radius: 6px;
    padding: 5px 8px;
    position: absolute;
    z-index: 1000;
    bottom: 125%;
    left: 50%;
    transform: translateX(-50%);
    transition: opacity 0.2s;
    font-size: 0.9em;
    font-weight: normal;
    white-space: nowrap;
    line-height: 1.4;
    color: #fff !important;
}
.custom-tooltip::after {
    content: '';
    position: absolute;
    top: 100%;
    left: 50%;
    margin-left: -5px;
    border-width: 5px;
    border-style: solid;
    border-color: #333 transparent transparent transparent;
}
.aa-residue:hover .custom-tooltip {
    visibility: visible;
    opacity: 1;
}
.aa-hydrophobic { background-color: #dbeafe; color: #1e40af; }
.aa-polar { background-color: #dcfce7; color: #166534; }
.aa-acidic { background-color: #fee2e2; color: #991b1b; }
.aa-basic { background-color: #f3e8ff; color: #581c87; }
.legend-container {
    padding: 0px 10px 5px 10px;
}
.legend-item {
    display: inline-flex;
    align-items: center;
    margin-right: 15px;
    margin-bottom: 5px;
}
.legend-box {
    width: 15px;
    height: 15px;
    border: 1px solid #888;
    display: inline-block;
    margin-right: 5px;
    border-radius: 2px;
}
.legend-label {
    display: inline-block;
    font-size: 0.9em;
}
.sequence-container {
    border: 1px solid #ddd;
    border-radius: 8px;
    padding: 10px 50px;
    background-color: #fafafa;
    word-break: break-all;
    line-height: 1.2;
}
</style>
"""

# ===================================================================
# SECTION 2: DATA LOADING FUNCTIONS
# ===================================================================

@st.cache_data
def load_protein_lengths():
    """Load protein lengths CSV."""
    try:
        df = pd.read_csv("data/protein_lengths.csv")
        if 'accession' not in df.columns and 'raw_id' in df.columns:
            df['accession'] = df['raw_id'].str.extract(r'\|([^|]+)\|', expand=False)
        return df
    except Exception as e:
        st.error(f"Error loading protein_lengths.csv: {str(e)}")
        return pd.DataFrame()

@st.cache_data
def load_protein_annotations():
    """Load protein annotations CSV."""
    try:
        df = pd.read_csv("data/protein_annotation.csv")
        if 'ProteinID_acc' not in df.columns and 'ProteinID' in df.columns:
            df['ProteinID_acc'] = df['ProteinID'].str.extract(r'\|([^|]+)\|', expand=False)
        return df
    except Exception as e:
        st.error(f"Error loading protein_annotation.csv: {str(e)}")
        return pd.DataFrame()

@st.cache_data
def load_fasta_sequences():
    """Load FASTA sequences."""
    try:
        fasta_path = "data/protein_sequences.fasta"
        if not os.path.exists(fasta_path):
            return {}
        with open(fasta_path, 'r') as f:
            lines = f.readlines()
        return parse_fasta(lines)
    except Exception as e:
        st.warning(f"Could not load FASTA file: {str(e)}")
        return {}

def get_name_map(protein_lengths_df):
    """Create mapping from raw_id to display name."""
    if protein_lengths_df.empty:
        return {}
    
    if 'real_name' in protein_lengths_df.columns:
        labels = protein_lengths_df.apply(
            lambda row: row['real_name'] if pd.notna(row['real_name']) and row['real_name'] != "" else row['short_name'],
            axis=1
        )
    else:
        labels = protein_lengths_df['short_name']
    
    return dict(zip(protein_lengths_df['raw_id'], labels))

def load_dataset(dataset_key, protein_lengths_df, protein_annotations_df):
    """Load dataset based on the selected key."""
    if dataset_key not in FILE_INFO_LIST:
        return None, None
    
    info = FILE_INFO_LIST[dataset_key]
    stem = info['stem']
    dataset_type = info['type']
    
    try:
        if dataset_type == "standard":
            links_path = f"data/{stem}_Links.csv"
            if not os.path.exists(links_path):
                return None, None
            
            links = pd.read_csv(links_path)
            annots = protein_annotations_df.copy()
            
            # Extract accessions
            links['Protein1_acc'] = links['Protein1'].str.extract(r'\|([^|]+)\|', expand=False)
            links['Protein2_acc'] = links['Protein2'].str.extract(r'\|([^|]+)\|', expand=False)
            
        elif dataset_type == "standard_shortname":
            links_path = f"data/{stem}.csv"
            if not os.path.exists(links_path):
                return None, None
            
            links_orig = pd.read_csv(links_path)
            
            # Map short names to raw_ids
            name_to_rawid = protein_lengths_df[['short_name', 'raw_id']].drop_duplicates()
            
            links_mapped = links_orig.merge(
                name_to_rawid, left_on='Protein1', right_on='short_name', how='left'
            ).rename(columns={'raw_id': 'Protein1_raw'}).merge(
                name_to_rawid, left_on='Protein2', right_on='short_name', how='left'
            ).rename(columns={'raw_id': 'Protein2_raw'})
            
            links_mapped = links_mapped.dropna(subset=['Protein1_raw', 'Protein2_raw'])
            
            # Drop original Protein1 and Protein2 columns, rename _raw columns
            links = links_mapped.drop(columns=['Protein1', 'Protein2']).copy()
            links = links.rename(columns={'Protein1_raw': 'Protein1', 'Protein2_raw': 'Protein2'})
            
            links['Protein1_acc'] = links['Protein1'].str.extract(r'\|([^|]+)\|', expand=False)
            links['Protein2_acc'] = links['Protein2'].str.extract(r'\|([^|]+)\|', expand=False)
            
            annots = protein_annotations_df.copy()
            
        else:
            # Proxl format (not fully implemented in original for these datasets)
            return None, None
        
        return links, annots
    
    except Exception as e:
        st.error(f"Error loading dataset {dataset_key}: {str(e)}")
        return None, None

# ===================================================================
# SECTION 3: VISUALIZATION FUNCTIONS
# ===================================================================

def generate_linear_plot(plot_data, show_intra_links=True):
    """Generate publication-quality interactive linear plot using Plotly."""
    if plot_data is None or plot_data['sector_df'].empty:
        return None
    
    sector_df = plot_data['sector_df'].copy()
    links_df = plot_data['links_df'].copy()
    annotations_bed = plot_data['annotations_bed'].copy()
    color_palette = plot_data['color_palette']
    
    # Sort proteins by length (descending)
    sector_df = sector_df.sort_values('end', ascending=False)
    protein_order = sector_df['name'].tolist()
    
    # Create figure with white template
    fig = go.Figure()
    
    max_len = sector_df['end'].max()
    
    # Create numeric y positions with increased spacing
    y_positions = {name: idx * 2 for idx, name in enumerate(protein_order)}  # Increased spacing (2x)
    
    # Add protein tracks as thin rectangles using shapes
    for idx, (_, row) in enumerate(sector_df.iterrows()):
        protein_name = row['name']
        color = color_palette.get(protein_name, '#808080')
        y_pos = y_positions[protein_name]
        
        # Use shapes for thinner, more precise rectangles
        fig.add_shape(
            type="rect",
            x0=row['start'],
            x1=row['end'],
            y0=y_pos - 0.2,  # Thinner track
            y1=y_pos + 0.2,
            fillcolor=color,
            line=dict(color='black', width=1.5),
            layer='below'
        )
        
        # Add protein name label
        fig.add_annotation(
            x=max_len * 0.01,
            y=y_pos,
            text=f"<b>{protein_name}</b>",
            showarrow=False,
            font=dict(size=11, color='black'),
            xanchor='left',
            bgcolor='white',
            bordercolor='black',
            borderwidth=1
        )
    
    # Add annotations (domain blocks) as overlays
    if not annotations_bed.empty:
        for _, row in annotations_bed.iterrows():
            protein_name = row['chr']
            if protein_name in y_positions:
                y_pos = y_positions[protein_name]
                color = color_palette.get(protein_name, '#808080')
                fig.add_shape(
                    type="rect",
                    x0=row['start'],
                    x1=row['end'],
                    y0=y_pos - 0.2,
                    y1=y_pos + 0.2,
                    fillcolor=color,
                    opacity=0.6,
                    line=dict(color='#666666', width=1),
                    layer='below'
                )
    
    # Process links
    if 'NumPSMs' in links_df.columns:
        links_df['AlphaValue'] = links_df['NumPSMs']
        score_col = 'NumPSMs'
    elif 'Score' in links_df.columns:
        links_df['AlphaValue'] = links_df['Score']
        score_col = 'Score'
    else:
        links_df['AlphaValue'] = 1.0
        score_col = None
    
    # Normalize alpha values for better visibility
    if links_df['AlphaValue'].max() > links_df['AlphaValue'].min():
        min_val = links_df['AlphaValue'].min()
        max_val = links_df['AlphaValue'].max()
        links_df['AlphaValue'] = 0.3 + 0.7 * ((links_df['AlphaValue'] - min_val) / (max_val - min_val))
    else:
        links_df['AlphaValue'] = 0.7
    
    # Separate inter and intra links
    inter_links = links_df[links_df['P1_clean'] != links_df['P2_clean']].copy()
    intra_links = links_df[
        (links_df['P1_clean'] == links_df['P2_clean']) & 
        (links_df['LinkPos1'] != links_df['LinkPos2'])
    ].copy()
    
    # Add inter-links (thin, transparent straight lines)
    for _, link in inter_links.iterrows():
        p1 = link['P1_clean']
        p2 = link['P2_clean']
        
        if p1 in y_positions and p2 in y_positions:
            color = color_palette.get(p1, '#808080')
            # Force opacity to 0.3 for inter-links
            alpha = 0.3
            score_text = f" | Score: {link[score_col]:.2f}" if score_col else ""
            
            fig.add_trace(go.Scatter(
                x=[link['LinkPos1'], link['LinkPos2']],
                y=[y_positions[p1] + 0.2, y_positions[p2] + 0.2],
                mode='lines',
                line=dict(color=color, width=1),  # Thinner lines
                opacity=alpha,
                showlegend=False,
                hovertemplate=f"<b>{p1}</b> (Pos {link['LinkPos1']}) <-> <b>{p2}</b> (Pos {link['LinkPos2']}){score_text}<extra></extra>"
            ))
    
    # Add intra-links (Bézier curves) if enabled
    if show_intra_links and not intra_links.empty:
        for _, link in intra_links.iterrows():
            p1 = link['P1_clean']
            
            if p1 in y_positions:
                y_pos = y_positions[p1]
                color = color_palette.get(p1, '#808080')
                alpha = link['AlphaValue']
                score_text = f" | Score: {link[score_col]:.2f}" if score_col else ""
                
                # Create Bézier curve for smooth arc (loops upward)
                x_start = link['LinkPos1']
                x_end = link['LinkPos2']
                
                # Control points for cubic Bézier curve
                control_offset = abs(x_end - x_start) * 0.4  # Horizontal offset for control points
                arc_height = abs(x_end - x_start) * 0.2  # Height of arc
                
                # Generate Bézier curve points
                n_points = 50
                t = np.linspace(0, 1, n_points)
                
                # Cubic Bézier: P(t) = (1-t)³P₀ + 3(1-t)²tP₁ + 3(1-t)t²P₂ + t³P₃
                x_curve = ((1-t)**3 * x_start + 
                          3*(1-t)**2*t * (x_start + control_offset) + 
                          3*(1-t)*t**2 * (x_end - control_offset) + 
                          t**3 * x_end)
                
                # Y curve loops upward
                y_curve = ((1-t)**3 * (y_pos + 0.2) + 
                          3*(1-t)**2*t * (y_pos + 0.2 + arc_height) + 
                          3*(1-t)*t**2 * (y_pos + 0.2 + arc_height) + 
                          t**3 * (y_pos + 0.2))
                
                fig.add_trace(go.Scatter(
                    x=x_curve,
                    y=y_curve,
                    mode='lines',
                    line=dict(color=color, width=2, smoothing=1.3),
                    opacity=alpha,
                    showlegend=False,
                    hovertemplate=f"<b>{p1}</b> (Pos {link['LinkPos1']} <-> Pos {link['LinkPos2']}){score_text}<extra></extra>"
                ))
    
    # Update layout for publication quality
    fig.update_layout(
        template="plotly_white",
        height=500,
        xaxis=dict(
            title=dict(text="<b>Amino Acid Position</b>", font=dict(size=14)),
            range=[-max_len * 0.02, max_len * 1.02],
            showgrid=False,
            showline=True,
            linecolor="black",
            linewidth=2,
            mirror=True,
            tickfont=dict(size=12)
        ),
        yaxis=dict(
            title="",
            tickmode='array',
            tickvals=list(y_positions.values()),
            ticktext=[f"<b>{p}</b>" for p in protein_order],
            showgrid=False,
            showline=False,
            showticklabels=False,  # Hide default labels since we use annotations
            tickfont=dict(size=12, color='black'),
            automargin=True
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin=dict(l=30, r=20, t=50, b=50),
        hovermode='closest'
    )
    
    return fig

def generate_circos_plot(plot_data, show_intra_links=True):
    """Generate interactive circos plot using Plotly polar coordinates."""
    if plot_data is None or plot_data['sector_df'].empty:
        return None
    
    sector_df = plot_data['sector_df'].copy()
    links_df = plot_data['links_df'].copy()
    annotations_bed = plot_data['annotations_bed'].copy()
    color_palette = plot_data['color_palette']
    
    # Calculate total length for angle distribution
    total_length = sector_df['end'].sum()
    gap_degrees = 0.5  # Gap between sectors in degrees
    
    # Create angle mapping for each protein
    current_angle = 90  # Start at top (90 degrees)
    sector_angles = {}
    sector_ranges = {}
    
    for _, row in sector_df.iterrows():
        protein_name = row['name']
        length = row['end']
        
        # Calculate angle span for this protein
        angle_span = (length / total_length) * (360 - gap_degrees * len(sector_df))
        
        # Store angle range
        sector_ranges[protein_name] = {
            'start_angle': current_angle,
            'end_angle': current_angle + angle_span,
            'length': length
        }
        
        # Move to next sector
        current_angle += angle_span + gap_degrees
    
    # Create position to angle mapping functions
    for protein_name, range_info in sector_ranges.items():
        def make_pos_to_angle(protein, start_angle, end_angle, length):
            def pos_to_angle(pos):
                angle_range = end_angle - start_angle
                normalized_pos = pos / length
                return start_angle + normalized_pos * angle_range
            return pos_to_angle
        
        sector_angles[protein_name] = make_pos_to_angle(
            protein_name,
            range_info['start_angle'],
            range_info['end_angle'],
            range_info['length']
        )
    
    # Create figure with polar subplot
    fig = go.Figure()
    
    # Define radii for different tracks
    outer_radius = 1.0  # Outer track for protein names
    inner_radius = 0.85  # Inner track for annotations
    link_radius = 0.75  # Radius for links
    
    # Add protein sectors (outer track)
    for protein_name, range_info in sector_ranges.items():
        color = color_palette.get(protein_name, '#808080')
        
        # Create arc for protein sector
        theta_start = range_info['start_angle']
        theta_end = range_info['end_angle']
        theta_range_deg = np.linspace(theta_start, theta_end, 100)
        
        # Outer track (protein name background) - create filled sector
        # Need to close the sector by adding center point
        theta_with_center = np.concatenate([[theta_start], theta_range_deg, [theta_end], [theta_start]])
        r_with_center = np.concatenate([[0.01], [outer_radius] * 100, [0.01], [0.01]])
        
        fig.add_trace(go.Scatterpolar(
            r=r_with_center,
            theta=theta_with_center,
            mode='lines',
            line=dict(color=color, width=2),
            fill='toself',
            fillcolor=color,
            showlegend=False,
            hoverinfo='skip',
            name=protein_name
        ))
        
        # Add protein name label at middle of sector
        mid_angle_deg = (theta_start + theta_end) / 2
        label_radius = outer_radius * 0.5  # Place label at middle of sector
        
        # Add text annotation using Scatterpolar with text
        fig.add_trace(go.Scatterpolar(
            r=[label_radius],
            theta=[mid_angle_deg],
            mode='text',
            text=[protein_name],
            textfont=dict(size=11, color='white', family='Arial Black'),
            showlegend=False,
            hoverinfo='skip'
        ))
    
    # Add annotation tracks
    if not annotations_bed.empty:
        for _, annot in annotations_bed.iterrows():
            protein_name = annot['chr']
            if protein_name not in sector_ranges:
                continue
            
            range_info = sector_ranges[protein_name]
            pos_to_angle = sector_angles[protein_name]
            
            start_angle = pos_to_angle(annot['start'])
            end_angle = pos_to_angle(annot['end'])
            
            color = color_palette.get(protein_name, '#808080')
            
            # Create annotation arc
            theta_range_deg = np.linspace(start_angle, end_angle, 50)
            
            # Close the arc to create a filled shape
            theta_closed = np.concatenate([[start_angle], theta_range_deg, [end_angle], [start_angle]])
            r_closed = np.concatenate([[inner_radius - 0.02], [inner_radius] * 50, [inner_radius - 0.02], [inner_radius - 0.02]])
            
            fig.add_trace(go.Scatterpolar(
                r=r_closed,
                theta=theta_closed,
                mode='lines',
                line=dict(color=color, width=15),
                fill='toself',
                fillcolor=color,
                opacity=0.6,
                showlegend=False,
                hovertemplate=f"{protein_name}: {annot['start']}-{annot['end']} ({annot['name']})<extra></extra>"
            ))
    
    # Process links
    if 'NumPSMs' in links_df.columns:
        links_df['AlphaValue'] = links_df['NumPSMs']
    elif 'Score' in links_df.columns:
        links_df['AlphaValue'] = links_df['Score']
    else:
        links_df['AlphaValue'] = 1.0
    
    # Normalize alpha values
    if links_df['AlphaValue'].max() > links_df['AlphaValue'].min():
        min_val = links_df['AlphaValue'].min()
        max_val = links_df['AlphaValue'].max()
        links_df['AlphaValue'] = 0.1 + 0.9 * ((links_df['AlphaValue'] - min_val) / (max_val - min_val))
    else:
        links_df['AlphaValue'] = 0.8
    
    # Filter links
    links_to_plot = links_df.copy()
    if not show_intra_links:
        links_to_plot = links_to_plot[links_to_plot['P1_clean'] != links_to_plot['P2_clean']]
    
    # Add links
    for _, link in links_to_plot.iterrows():
        p1 = link['P1_clean']
        p2 = link['P2_clean']
        pos1 = link['LinkPos1']
        pos2 = link['LinkPos2']
        
        if p1 not in sector_angles or p2 not in sector_angles:
            continue
        
        # Get angles for both positions
        angle1 = sector_angles[p1](pos1)
        angle2 = sector_angles[p2](pos2)
        
        # Convert to radians
        theta1 = np.radians(angle1)
        theta2 = np.radians(angle2)
        
        # Get color
        link_color = color_palette.get(p1, '#808080')
        alpha = link['AlphaValue']
        
        # Create curved path for link (bezier-like curve)
        # Use multiple points to create a smooth arc
        n_points = 50
        theta_points = np.linspace(theta1, theta2, n_points)
        
        # Create curved radius (starts and ends at link_radius, curves outward)
        curve_height = 0.15  # How much the curve extends outward
        t = np.linspace(0, 1, n_points)
        # Quadratic curve: r = link_radius + curve_height * 4 * t * (1 - t)
        r_points = link_radius + curve_height * 4 * t * (1 - t)
        
        # Convert rgba color with alpha
        def hex_to_rgba(hex_color, alpha):
            hex_color = hex_color.lstrip('#')
            r = int(hex_color[0:2], 16)
            g = int(hex_color[2:4], 16)
            b = int(hex_color[4:6], 16)
            return f'rgba({r},{g},{b},{alpha})'
        
        rgba_color = hex_to_rgba(link_color, alpha)
        
        fig.add_trace(go.Scatterpolar(
            r=r_points,
            theta=np.degrees(theta_points),
            mode='lines',
            line=dict(color=rgba_color, width=2),
            showlegend=False,
            hovertemplate=f"{p1}({pos1}) - {p2}({pos2})<extra></extra>"
        ))
    
    # Add axis track (scale)
    # Add tick marks every 100 amino acids
    for protein_name, range_info in sector_ranges.items():
        length = range_info['length']
        pos_to_angle = sector_angles[protein_name]
        
        # Add ticks every 100 positions
        tick_positions = np.arange(0, length + 100, 100)
        for tick_pos in tick_positions:
            if tick_pos <= length:
                tick_angle = pos_to_angle(tick_pos)
                theta = np.radians(tick_angle)
                
                # Draw tick line
                fig.add_trace(go.Scatterpolar(
                    r=[outer_radius, outer_radius + 0.02],
                    theta=[tick_angle, tick_angle],
                    mode='lines',
                    line=dict(color='black', width=1),
                    showlegend=False,
                    hoverinfo='skip'
                ))
    
    # Update layout for polar plot
    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 1.1],
                showticklabels=False,
                showgrid=True,
                gridcolor='lightgray',
                gridwidth=1
            ),
            angularaxis=dict(
                rotation=90,  # Start at top
                direction='counterclockwise',
                showticklabels=False,
                showgrid=True,
                gridcolor='lightgray',
                gridwidth=1
            ),
            bgcolor='white'
        ),
        height=700,
        showlegend=False,
        plot_bgcolor='white'
    )
    
    return fig

def generate_lollipop_plot(plot_data, height_factor=0.6):
    """Generate lollipop plot using Plotly."""
    if plot_data is None or plot_data['plot_data'].empty:
        return None
    
    plot_df = plot_data['plot_data'].copy()
    target_length = plot_data['target_length']
    target_annotations = plot_data['target_annotations'].copy()
    color_palette = plot_data['color_palette']
    target_name = plot_data.get('target_name', 'Target Protein')
    
    fig = go.Figure()
    
    # Add protein backbone (grey rectangle)
    fig.add_shape(
        type="rect",
        x0=1, y0=0,
        x1=target_length, y1=1.5,
        fillcolor="grey",
        line=dict(color="black", width=1),
        opacity=0.8,
        layer="below"
    )
    
    # Add annotation blocks
    if not target_annotations.empty:
        for _, annot in target_annotations.iterrows():
            fig.add_shape(
                type="rect",
                x0=annot['StartRes'], y0=0,
                x1=annot['EndRes'], y1=1.5,
                fillcolor="#696969",
                line=dict(color="black", width=1),
                layer="below"
            )
            
            # Add annotation label
            mid_x = (annot['StartRes'] + annot['EndRes']) / 2
            fig.add_annotation(
                x=mid_x,
                y=0.75,
                text=annot['AnnotName'],
                showarrow=False,
                font=dict(size=10, color="white"),
                bgcolor="rgba(0,0,0,0)"
            )
    
    # Group by target position for lollipop stems
    for target_pos, group in plot_df.groupby('target_pos'):
        y_positions = group['y_position'].tolist()
        partner_names = group['partner_name'].tolist()
        partner_positions = group['partner_pos'].tolist()
        
        # Add stems (segments from y=0 to y_position)
        for y_pos, partner_name, partner_pos in zip(y_positions, partner_names, partner_positions):
            color = color_palette.get(partner_name, '#808080')
            
            # Stem
            fig.add_trace(go.Scatter(
                x=[target_pos, target_pos],
                y=[0, y_pos],
                mode='lines',
                line=dict(color=color, width=0.7),
                showlegend=False,
                hovertemplate=f"Target Position: {target_pos}<br>Partner: {partner_name} (Pos: {partner_pos})<extra></extra>"
            ))
            
            # Head (dot)
            fig.add_trace(go.Scatter(
                x=[target_pos],
                y=[y_pos],
                mode='markers',
                marker=dict(
                    color=color,
                    size=8,
                    opacity=0.8
                ),
                name=partner_name,
                showlegend=False,
                hovertemplate=f"Target Position: {target_pos}<br>Partner: {partner_name} (Pos: {partner_pos})<extra></extra>"
            ))
    
    # Calculate y-axis breaks
    counts = plot_df.groupby('target_pos').size()
    max_count = counts.max() if not counts.empty else 0
    y_axis_breaks = list(range(0, max_count + 1, 5 if max_count > 50 else (2 if max_count > 20 else 1)))
    
    # Update layout
    fig.update_layout(
        title=dict(text=target_name, x=0.5, font=dict(size=20)),
        height=600,
        xaxis=dict(
            title="Amino Acid Position",
            range=[0, target_length + 100],
            showgrid=False,
            showline=True,
            linecolor="black"
        ),
        yaxis=dict(
            title="# Cross-links",
            tickmode='array',
            tickvals=[-b * height_factor for b in y_axis_breaks],
            ticktext=[str(b) for b in y_axis_breaks],
            range=[min(plot_df['y_position'].min() - 2 * height_factor, -0.5), 2.5],
            showgrid=False,
            showline=True,
            linecolor="black"
        ),
        plot_bgcolor='white',
        margin=dict(l=80, r=50, t=80, b=50)
    )
    
    return fig

# ===================================================================
# SECTION 4: MAIN APPLICATION
# ===================================================================

# Initialize session state
if 'plot_data_circos' not in st.session_state:
    st.session_state.plot_data_circos = None
if 'highlight_range' not in st.session_state:
    st.session_state.highlight_range = None

# Load master data
protein_lengths_master = load_protein_lengths()
master_annotations = load_protein_annotations()
loaded_sequences = load_fasta_sequences()

# Title
st.title("📊 2D Data Explorer")
st.markdown("Visualize cross-linking mass spectrometry data in 2D views")

# Sidebar
with st.sidebar:
    st.header("📁 Dataset Selection")
    
    dataset_options = {
        "Both CBAF and PBAF": "both_cbaf_pbaf",
        "Canonical (Between Subunits)": "canonical_between",
        "Canonical (Intra- and Inter-)": "canonical_intra",
        "cBAF + Nucleosome": "nucleosome",
        "PBAF Only": "pbaf",
        "ncBAF (Inter-subunit)": "ncbaf_inter1",
        "ncBAF (Intra-subunit)": "ncbaf_intra1",
        "cBAF GR (Proxl)": "cbaf_gr_proxl",
        "cBAF GR NCP (Proxl)": "cbaf_gr_ncp_proxl",
        "HA DPF2 (Proxl)": "ha_dpf2_proxl",
        "HA DPF2 NCP (Proxl)": "ha_dpf2_ncp_proxl",
        "SS18 Flag (Proxl)": "ss18_flag_proxl",
        "HA PHF10 (Proxl)": "phf10_ha_proxl",
        "HA BRD7 (Proxl)": "ha_brd7_proxl"
    }
    
    selected_dataset_display = st.selectbox(
        "1. Select Dataset",
        options=list(dataset_options.keys()),
        key="dataset_selector"
    )
    
    selected_dataset_key = dataset_options[selected_dataset_display]
    
    st.divider()
    
    st.header("🎨 Plot Controls")
    
    # Load dataset
    links_df_orig, annotations_df_orig = load_dataset(
        selected_dataset_key, protein_lengths_master, master_annotations
    )
    
    if links_df_orig is None:
        st.error("Could not load dataset. Please check that data files exist in the 'data/' directory.")
        st.stop()
    
    # Get available proteins
    all_protein_ids = set(links_df_orig['Protein1'].unique()) | set(links_df_orig['Protein2'].unique())
    if not annotations_df_orig.empty:
        all_protein_ids.update(annotations_df_orig['ProteinID'].unique())
    
    relevant_proteins = protein_lengths_master[
        protein_lengths_master['raw_id'].isin(all_protein_ids)
    ].copy()
    
    # Create name mapping
    name_map_dict = get_name_map(relevant_proteins)
    
    # Create selector options
    if 'real_name' in relevant_proteins.columns:
        display_names = relevant_proteins.apply(
            lambda row: row['real_name'] if pd.notna(row['real_name']) and row['real_name'] != "" else row['short_name'],
            axis=1
        )
    else:
        display_names = relevant_proteins['short_name']
    
    protein_options = dict(zip(display_names, relevant_proteins['raw_id']))
    
    selected_proteins = st.multiselect(
        "Select Proteins:",
        options=list(protein_options.keys()),
        default=list(protein_options.keys())[:min(10, len(protein_options))] if protein_options else [],
        key="selected_proteins_circos"
    )
    
    show_intra_links = st.checkbox(
        "Show Self-Links (Intra-Protein)",
        value=True,
        key="show_intra_links"
    )
    
    plot_button = st.button("Generate Plot", type="primary", use_container_width=True)
    
    st.divider()
    
    st.header("💾 Download Settings")
    download_width = st.number_input("Width (in):", min_value=1, max_value=50, value=10, key="download_width")
    download_height = st.number_input("Height (in):", min_value=1, max_value=50, value=10, key="download_height")

# Main content area
if plot_button or st.session_state.plot_data_circos is not None:
    if not selected_proteins:
        st.warning("Please select at least one protein.")
    else:
        # Process data for plotting
        selected_raw_ids = [protein_options[p] for p in selected_proteins]
        
        # Filter links
        links_df = links_df_orig[
            (links_df_orig['Protein1'].isin(selected_raw_ids)) &
            (links_df_orig['Protein2'].isin(selected_raw_ids))
        ].copy()
        
        # Filter annotations
        annotations_df = annotations_df_orig[
            annotations_df_orig['ProteinID'].isin(selected_raw_ids)
        ].copy() if not annotations_df_orig.empty else pd.DataFrame()
        
        # Get protein lengths
        protein_lengths = protein_lengths_master[
            protein_lengths_master['raw_id'].isin(selected_raw_ids)
        ].copy()
        
        if protein_lengths.empty:
            st.error("No valid proteins selected.")
            st.stop()
        
        # Create clean names
        protein_lengths['clean_name'] = protein_lengths['raw_id'].map(name_map_dict)
        protein_lengths = protein_lengths.dropna(subset=['clean_name'])
        
        # Create sector dataframe
        sector_df = protein_lengths[['clean_name', 'length']].copy()
        sector_df.columns = ['name', 'end']
        sector_df['start'] = 0
        
        # Map links to clean names
        id_map = dict(zip(protein_lengths['raw_id'], protein_lengths['clean_name']))
        links_df['P1_clean'] = links_df['Protein1'].map(id_map)
        links_df['P2_clean'] = links_df['Protein2'].map(id_map)
        links_df = links_df.dropna(subset=['P1_clean', 'P2_clean'])
        
        # Prepare annotations
        annotations_df['clean_name'] = annotations_df['ProteinID'].map(id_map)
        annotations_bed = annotations_df.dropna(subset=['clean_name'])[
            ['clean_name', 'StartRes', 'EndRes', 'AnnotName']
        ].copy()
        annotations_bed.columns = ['chr', 'start', 'end', 'name']
        
        # Create color palette
        all_proteins_in_plot = sorted(sector_df['name'].unique())
        color_palette = {}
        for protein in all_proteins_in_plot:
            if protein in FIXED_PALETTE:
                color_palette[protein] = FIXED_PALETTE[protein]
            else:
                # Use a default color scheme for unmapped proteins
                color_palette[protein] = px.colors.qualitative.Set3[
                    len([p for p in all_proteins_in_plot if p in FIXED_PALETTE]) % len(px.colors.qualitative.Set3)
                ]
        
        # Store plot data
        plot_data = {
            'links_df': links_df,
            'sector_df': sector_df,
            'annotations_bed': annotations_bed,
            'annotations_df': annotations_df,
            'color_palette': color_palette
        }
        
        if plot_button:
            st.session_state.plot_data_circos = plot_data
        
        # Create tabs for different views
        tab1, tab2 = st.tabs(["Global View", "Lollipop View"])
        
        with tab1:
            st.header("Global View")
            
            if st.session_state.plot_data_circos:
                # Generate Linear plot once (reused for display and download)
                fig_linear = generate_linear_plot(
                    st.session_state.plot_data_circos,
                    show_intra_links=show_intra_links
                )
                
                # Create side-by-side columns: Circos (left) and Linear (right)
                col1, col2 = st.columns([1, 1])
                
                with col1:
                    st.subheader("Circos Plot")
                    # Generate Circos plot
                    plot_data = st.session_state.plot_data_circos
                    links_df = plot_data['links_df'].copy()
                    sector_df = plot_data['sector_df'].copy()
                    annotations_bed = plot_data['annotations_bed'].copy()
                    
                    # Filter intra-links if needed
                    if not show_intra_links:
                        links_df = links_df[links_df['P1_clean'] != links_df['P2_clean']].copy()
                    
                    # Ensure temp directory exists
                    temp_dir = "temp"
                    os.makedirs(temp_dir, exist_ok=True)
                    
                    # Prepare file paths
                    links_file = os.path.join(temp_dir, "temp_links.csv")
                    sectors_file = os.path.join(temp_dir, "temp_sectors.csv")
                    annots_file = os.path.join(temp_dir, "temp_annots.csv")
                    output_file = os.path.join(temp_dir, "circos_plot.png")
                    
                    # Save dataframes to CSV
                    try:
                        # Prepare links_df for R (ensure required columns exist)
                        links_export = links_df[['P1_clean', 'LinkPos1', 'P2_clean', 'LinkPos2']].copy()
                        if 'NumPSMs' in links_df.columns:
                            links_export['NumPSMs'] = links_df['NumPSMs']
                        elif 'Score' in links_df.columns:
                            links_export['Score'] = links_df['Score']
                        else:
                            links_export['Score'] = 1.0
                        
                        links_export.to_csv(links_file, index=False)
                        
                        # Prepare sector_df for R
                        sector_export = sector_df[['name', 'start', 'end']].copy()
                        sector_export.to_csv(sectors_file, index=False)
                        
                        # Prepare annotations_bed for R (rename columns to match R expectations)
                        if not annotations_bed.empty:
                            annots_export = annotations_bed[['chr', 'start', 'end', 'name']].copy()
                            annots_export.to_csv(annots_file, index=False)
                        else:
                            # Create empty file with correct structure
                            pd.DataFrame(columns=['chr', 'start', 'end', 'name']).to_csv(annots_file, index=False)
                        
                        # Call R script
                        r_script_path = "r_scripts/generate_circos.R"
                        
                        if not os.path.exists(r_script_path):
                            st.error(f"R script not found at {r_script_path}. Please ensure the script exists.")
                        else:
                            # Run R script
                            cmd = [
                                "Rscript",
                                r_script_path,
                                links_file,
                                sectors_file,
                                annots_file,
                                output_file
                            ]
                            
                            with st.spinner("Generating Circos plot..."):
                                result = subprocess.run(
                                    cmd,
                                    capture_output=True,
                                    text=True,
                                    timeout=60  # 60 second timeout
                                )
                            
                            if result.returncode == 0:
                                # Check if output file was created
                                if os.path.exists(output_file):
                                    st.image(output_file, use_container_width=True)
                                    
                                    # Download buttons for Circos Plot
                                    with open(output_file, 'rb') as f:
                                        circos_img_bytes = f.read()
                                    
                                    circos_dl_col1, circos_dl_col2 = st.columns(2)
                                    with circos_dl_col1:
                                        st.download_button(
                                            label="📥 Download as PNG",
                                            data=circos_img_bytes,
                                            file_name=f"circos_plot_{selected_dataset_key}.png",
                                            mime="image/png"
                                        )
                                    with circos_dl_col2:
                                        # Convert PNG to PDF using PIL if available
                                        try:
                                            from PIL import Image
                                            import io
                                            img = Image.open(io.BytesIO(circos_img_bytes))
                                            pdf_buffer = io.BytesIO()
                                            # Convert to RGB if necessary (PDF doesn't support RGBA)
                                            if img.mode == 'RGBA':
                                                rgb_img = Image.new('RGB', img.size, (255, 255, 255))
                                                rgb_img.paste(img, mask=img.split()[3])
                                                img = rgb_img
                                            img.save(pdf_buffer, format='PDF')
                                            pdf_bytes = pdf_buffer.getvalue()
                                            st.download_button(
                                                label="📥 Download as PDF",
                                                data=pdf_bytes,
                                                file_name=f"circos_plot_{selected_dataset_key}.pdf",
                                                mime="application/pdf"
                                            )
                                        except ImportError:
                                            st.info("PDF requires PIL/Pillow")
                                        except Exception:
                                            st.info("PDF conversion failed")
                                else:
                                    st.error(f"R script completed but output file not found at {output_file}")
                                    if result.stderr:
                                        st.code(result.stderr, language='text')
                            else:
                                st.error("Error generating Circos plot.")
                                if result.stderr:
                                    st.code(result.stderr, language='text')
                    
                    except Exception as e:
                        st.error(f"Error: {str(e)}")
                        import traceback
                        st.code(traceback.format_exc(), language='text')
                
                with col2:
                    st.subheader("Interactive Linear Plot")
                    if fig_linear:
                        st.plotly_chart(fig_linear, use_container_width=True)
                        
                        # Download buttons for Linear Plot
                        linear_dl_col1, linear_dl_col2 = st.columns(2)
                        with linear_dl_col1:
                            try:
                                img_bytes = fig_linear.to_image(format="png", width=1200, height=800)
                                st.download_button(
                                    label="📥 Download as PNG",
                                    data=img_bytes,
                                    file_name=f"linear_plot_{selected_dataset_key}.png",
                                    mime="image/png"
                                )
                            except Exception:
                                st.info("PNG requires kaleido")
                        with linear_dl_col2:
                            try:
                                pdf_bytes = fig_linear.to_image(format="pdf", width=1200, height=800)
                                st.download_button(
                                    label="📥 Download as PDF",
                                    data=pdf_bytes,
                                    file_name=f"linear_plot_{selected_dataset_key}.pdf",
                                    mime="application/pdf"
                                )
                            except Exception:
                                st.info("PDF requires kaleido")
                    else:
                        st.warning("Could not generate linear plot.")
            else:
                st.info("Click 'Generate Plot' to create visualizations.")
            
            st.divider()
            
            st.subheader("Crosslink Details")
            if st.session_state.plot_data_circos and not st.session_state.plot_data_circos['links_df'].empty:
                display_df = st.session_state.plot_data_circos['links_df'][
                    ['P1_clean', 'LinkPos1', 'P2_clean', 'LinkPos2']
                ].copy()
                display_df.columns = ['Protein 1', 'Position 1', 'Protein 2', 'Position 2']
                
                if 'NumPSMs' in st.session_state.plot_data_circos['links_df'].columns:
                    display_df['PSMs'] = st.session_state.plot_data_circos['links_df']['NumPSMs']
                elif 'Score' in st.session_state.plot_data_circos['links_df'].columns:
                    display_df['Score'] = st.session_state.plot_data_circos['links_df']['Score']
                
                st.dataframe(display_df, width='stretch', hide_index=True)
            else:
                st.info("No crosslinks to display.")
            
            st.divider()
            
            st.subheader("Annotation Details")
            if st.session_state.plot_data_circos and not st.session_state.plot_data_circos['annotations_df'].empty:
                annot_display = st.session_state.plot_data_circos['annotations_df'][
                    ['clean_name', 'AnnotName', 'StartRes', 'EndRes']
                ].copy()
                annot_display.columns = ['Protein', 'Annotation Name', 'Start', 'End']
                st.dataframe(annot_display, width='stretch', hide_index=True)
            else:
                st.info("No annotations available.")
        
        with tab2:
            st.header("Lollipop View")
            
            # Check if dataset type is supported
            dataset_type = FILE_INFO_LIST[selected_dataset_key]['type']
            if dataset_type not in ['standard', 'standard_shortname']:
                st.warning("Lollipop View is not supported for Proxl datasets.")
            else:
                # Target protein selector
                target_protein_display = st.selectbox(
                    "Select Target Protein:",
                    options=list(protein_options.keys()),
                    key="target_protein_id_lollipop"
                )
                
                if target_protein_display:
                    target_protein_id = protein_options[target_protein_display]
                    target_display_name = name_map_dict.get(target_protein_id, target_protein_id)
                    
                    # Get target protein data
                    target_length = protein_lengths_master[
                        protein_lengths_master['raw_id'] == target_protein_id
                    ]['length'].values
                    
                    if len(target_length) > 0:
                        target_length = target_length[0]
                        
                        # Filter links for target protein
                        target_links = links_df_orig[
                            (links_df_orig['Protein1'] == target_protein_id) |
                            (links_df_orig['Protein2'] == target_protein_id)
                        ].copy()
                        
                        if not target_links.empty:
                            # Process links
                            processed_links = target_links.copy()
                            processed_links['target_pos'] = processed_links.apply(
                                lambda row: row['LinkPos1'] if row['Protein1'] == target_protein_id else row['LinkPos2'],
                                axis=1
                            )
                            processed_links['partner_id'] = processed_links.apply(
                                lambda row: row['Protein2'] if row['Protein1'] == target_protein_id else row['Protein1'],
                                axis=1
                            )
                            processed_links['partner_pos'] = processed_links.apply(
                                lambda row: row['LinkPos2'] if row['Protein1'] == target_protein_id else row['LinkPos1'],
                                axis=1
                            )
                            
                            # Get partner names
                            processed_links['partner_name'] = processed_links['partner_id'].map(name_map_dict)
                            # For self-links, use target name
                            processed_links.loc[
                                processed_links['partner_id'] == target_protein_id,
                                'partner_name'
                            ] = target_display_name
                            processed_links['partner_name'] = processed_links['partner_name'].fillna(
                                processed_links['partner_id'].str.extract(r'\|([^|]+)\|', expand=False)
                            )
                            
                            # Calculate y positions
                            height_factor = st.slider(
                                "Lollipop Height Factor:",
                                min_value=0.2,
                                max_value=2.0,
                                value=0.6,
                                step=0.1,
                                key="height_factor_lollipop"
                            )
                            
                            processed_links = processed_links.sort_values(['target_pos', 'partner_name'])
                            processed_links['y_position'] = processed_links.groupby('target_pos').cumcount() * -height_factor
                            
                            # Get target annotations
                            target_annotations = annotations_df_orig[
                                annotations_df_orig['ProteinID'] == target_protein_id
                            ].copy() if not annotations_df_orig.empty else pd.DataFrame()
                            
                            # Create color palette for partners
                            all_partners = sorted(processed_links['partner_name'].unique())
                            partner_color_palette = {}
                            for partner in all_partners:
                                if partner in FIXED_PALETTE:
                                    partner_color_palette[partner] = FIXED_PALETTE[partner]
                                else:
                                    partner_color_palette[partner] = px.colors.qualitative.Set3[
                                        len([p for p in all_partners if p in FIXED_PALETTE]) % len(px.colors.qualitative.Set3)
                                    ]
                            
                            # Prepare plot data
                            lollipop_plot_data = {
                                'plot_data': processed_links,
                                'target_length': target_length,
                                'target_annotations': target_annotations,
                                'color_palette': partner_color_palette,
                                'target_name': target_display_name
                            }
                            
                            # Generate plot
                            fig_lollipop = generate_lollipop_plot(lollipop_plot_data, height_factor=height_factor)
                            
                            if fig_lollipop:
                                st.plotly_chart(fig_lollipop, use_container_width=True)
                                
                                # Download buttons for Lollipop Plot
                                col_dl1, col_dl2 = st.columns(2)
                                with col_dl1:
                                    try:
                                        img_bytes = fig_lollipop.to_image(format="png", width=1200, height=800)
                                        st.download_button(
                                            label="📥 Download as PNG",
                                            data=img_bytes,
                                            file_name=f"lollipop_{target_display_name}_{selected_dataset_key}.png",
                                            mime="image/png"
                                        )
                                    except Exception:
                                        st.info("PNG download requires kaleido. Install with: pip install kaleido")
                                with col_dl2:
                                    try:
                                        pdf_bytes = fig_lollipop.to_image(format="pdf", width=1200, height=800)
                                        st.download_button(
                                            label="📥 Download as PDF",
                                            data=pdf_bytes,
                                            file_name=f"lollipop_{target_display_name}_{selected_dataset_key}.pdf",
                                            mime="application/pdf"
                                        )
                                    except Exception:
                                        st.info("PDF download requires kaleido. Install with: pip install kaleido")
                            
                            st.divider()
                            
                            # Sequence viewer
                            st.subheader("Interactive Protein Sequence")
                            
                            # Legend
                            st.markdown("""
                            <div class="legend-container">
                                <span class="legend-item">
                                    <span class="legend-box aa-hydrophobic"></span>
                                    <span class="legend-label">Hydrophobic</span>
                                </span>
                                <span class="legend-item">
                                    <span class="legend-box aa-polar"></span>
                                    <span class="legend-label">Polar</span>
                                </span>
                                <span class="legend-item">
                                    <span class="legend-box aa-acidic"></span>
                                    <span class="legend-label">Acidic</span>
                                </span>
                                <span class="legend-item">
                                    <span class="legend-box aa-basic"></span>
                                    <span class="legend-label">Basic</span>
                                </span>
                            </div>
                            """, unsafe_allow_html=True)
                            
                            # Highlight controls
                            col1, col2, col3 = st.columns([2, 2, 3])
                            with col1:
                                start_pos = st.number_input("Start:", min_value=1, max_value=target_length, value=1, key="start_pos")
                            with col2:
                                end_pos = st.number_input("End:", min_value=1, max_value=target_length, value=min(10, target_length), key="end_pos")
                            with col3:
                                if st.button("Highlight", key="highlight_btn"):
                                    st.session_state.highlight_range = (start_pos, end_pos)
                            
                            # Get sequence
                            target_accession = protein_lengths_master[
                                protein_lengths_master['raw_id'] == target_protein_id
                            ]['accession'].values
                            
                            if len(target_accession) > 0 and target_accession[0]:
                                target_accession = target_accession[0]
                                
                                # Find matching sequence
                                matching_header = None
                                for header in loaded_sequences.keys():
                                    if f"|{target_accession}|" in header:
                                        matching_header = header
                                        break
                                
                                if matching_header and matching_header in loaded_sequences:
                                    sequence = loaded_sequences[matching_header]
                                    
                                    # Generate HTML for sequence (no newlines)
                                    highlight_range = st.session_state.highlight_range
                                    
                                    spans = []
                                    for i, char in enumerate(sequence, 1):
                                        aa_class = get_aa_class(char)
                                        extra_class = ""
                                        if highlight_range and highlight_range[0] <= i <= highlight_range[1]:
                                            extra_class = "aa-highlight"
                                        
                                        title_text = f"{char} (Position: {i})"
                                        spans.append(f'<span class="aa-residue {aa_class} {extra_class}">{char}<span class="custom-tooltip">{title_text}</span></span>')
                                    
                                    seq_html = '<div class="sequence-container">' + ''.join(spans) + '</div>'
                                    
                                    st.markdown(SEQUENCE_VIEWER_CSS, unsafe_allow_html=True)
                                    st.markdown(seq_html, unsafe_allow_html=True)
                                else:
                                    st.info("Sequence not found in protein_sequences.fasta for the selected protein.")
                            else:
                                st.info("Could not find accession for selected protein.")
                            
                            st.divider()
                            
                            # Summary table - map partner_id to real_name
                            # Create mapping from raw_id to real_name
                            partner_id_to_real_name = {}
                            for _, row in protein_lengths_master.iterrows():
                                raw_id = row['raw_id']
                                if 'real_name' in row and pd.notna(row['real_name']) and row['real_name'] != "":
                                    partner_id_to_real_name[raw_id] = row['real_name']
                                elif 'short_name' in row:
                                    partner_id_to_real_name[raw_id] = row['short_name']
                            
                            # Map partner_id to real_name for summary
                            processed_links['partner_real_name'] = processed_links['partner_id'].map(partner_id_to_real_name)
                            # For self-links, use target real_name
                            target_real_name = partner_id_to_real_name.get(target_protein_id, target_display_name)
                            processed_links.loc[
                                processed_links['partner_id'] == target_protein_id,
                                'partner_real_name'
                            ] = target_real_name
                            # Fill any remaining NaN with partner_name
                            processed_links['partner_real_name'] = processed_links['partner_real_name'].fillna(processed_links['partner_name'])
                            
                            summary_df = processed_links.groupby('partner_real_name').size().reset_index(name='Number of Cross-links')
                            summary_df = summary_df.sort_values('Number of Cross-links', ascending=False)
                            summary_df.columns = ['Partner Subunit', 'Number of Cross-links']
                            st.dataframe(summary_df, width='stretch', hide_index=True)
                        else:
                            st.info(f"No cross-links found for {target_display_name} in the selected dataset.")
        

else:
    st.info("👈 Select a dataset and proteins in the sidebar, then click 'Generate Plot' to begin visualization.")

