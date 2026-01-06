import streamlit as st

# Page configuration
st.set_page_config(page_title="XL-MS Explorer", layout="wide")

# Custom CSS for modern card design
st.markdown("""
<style>
    .card {
        background-color: #ffffff;
        border-radius: 12px;
        padding: 1.5rem;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1), 0 1px 3px rgba(0, 0, 0, 0.08);
        transition: all 0.3s ease;
        border: 1px solid #e5e7eb;
        margin-bottom: 1.5rem;
        display: flex;
        flex-direction: column;
    }
    .card:hover {
        box-shadow: 0 10px 15px rgba(0, 0, 0, 0.15), 0 4px 6px rgba(0, 0, 0, 0.1);
        transform: translateY(-2px);
    }
    .card-title {
        font-size: 1.5rem;
        font-weight: 700;
        margin-bottom: 0.75rem;
        color: #1f2937;
    }
    .card-description {
        color: #4b5563;
        margin-bottom: 1rem;
        line-height: 1.6;
    }
    .card-features {
        margin-top: auto;
        padding-top: 1rem;
        border-top: 1px solid #e5e7eb;
    }
    .card-features-title {
        font-weight: 600;
        color: #374151;
        font-size: 0.9rem;
        margin-bottom: 0.5rem;
    }
    .card-features-list {
        color: #6b7280;
        font-size: 0.85rem;
        line-height: 1.8;
        margin: 0;
        padding-left: 1.2rem;
    }
    .gradient-text {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
        font-size: 3rem;
        font-weight: 800;
        text-align: center;
        margin-bottom: 1rem;
    }
    .subheader-text {
        text-align: left;
        color: #6b7280;
        font-size: 1.1rem;
        margin: 0 0 3rem 0;
        line-height: 1.7;
    }
</style>
""", unsafe_allow_html=True)

# Hero Section
st.markdown("""
<div class="gradient-text">Kadoch Lab Crosslink Explorer</div>
<p class="subheader-text">
This interactive dashboard is designed to visualize and analyze Cross-linking Mass Spectrometry (XL-MS) data for BAF complexes. 
It provides multiple perspectives on structural connectivity, allowing researchers to validate subunit interactions and map crosslinks to functional domains.
</p>
""", unsafe_allow_html=True)

# Feature Cards Stacked Vertically
st.markdown("""
<div class="card">
    <div class="card-title">🧬 3D Structure View (Structural Analysis)</div>
    <div class="card-description">
        Interact with 3D molecular structures (PDB) mapped with crosslinking data.
    </div>
    <div class="card-features">
        <div class="card-features-title">Features:</div>
        <ul class="card-features-list">
            <li><strong>Interactive Visualization:</strong> Rotate, zoom, and pan 3D protein structures using NGL.js viewer</li>
            <li><strong>Crosslink Mapping:</strong> Visualize crosslinks as colored cylinders connecting specific residue positions across protein chains</li>
            <li><strong>Filtering:</strong> Selectively view crosslinks based on confidence scores, distance thresholds, or specific protein-protein interactions</li>
            <li><strong>Link Inspector:</strong> Detailed table showing all visualized crosslinks with calculated distances (Å) and confidence scores</li>
        </ul>
    </div>
</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="card">
    <div class="card-title">📏 Linear View (Domain Mapping)</div>
    <div class="card-description">
        A 'genomic browser' style view that lays proteins out horizontally.
    </div>
    <div class="card-features">
        <div class="card-features-title">Features:</div>
        <ul class="card-features-list">
            <li><strong>Tracks:</strong> Each row represents a selected protein subunit, displayed as horizontal bars with distinct colors</li>
            <li><strong>Interactions:</strong> Arcs (Bézier curves) represent intra-protein links; straight lines represent inter-protein crosslinks</li>
            <li><strong>Interactive Hover:</strong> Hovering over any link reveals specific residue positions, protein names, and confidence scores</li>
            <li><strong>Domain Annotations:</strong> Functional domains and regions are overlaid as colored blocks on each protein track</li>
        </ul>
    </div>
</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="card">
    <div class="card-title">🍭 Lollipop View (Deep Dive)</div>
    <div class="card-description">
        Focuses on a single 'Target Protein' to visualize crosslink density.
    </div>
    <div class="card-features">
        <div class="card-features-title">Features:</div>
        <ul class="card-features-list">
            <li><strong>Stick Density:</strong> Each stick (lollipop) represents a crosslinked residue position, with height indicating the number of links at that position</li>
            <li><strong>Color Coding:</strong> The dot color at the top of each stick indicates the Partner Protein linked to this residue</li>
            <li><strong>Sequence Viewer:</strong> Full amino acid sequence displayed below the plot, chemically color-coded by residue type (Hydrophobic, Polar, Acidic, Basic)</li>
            <li><strong>Summary Table:</strong> Ranks all interacting partner proteins by crosslink count, using real_name identifiers (e.g., SMARCA4)</li>
            <li><strong>Interactive Highlighting:</strong> Select position ranges to highlight specific regions in the sequence viewer</li>
        </ul>
    </div>
</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="card">
    <div class="card-title">⭕ Circos View (Global Topology)</div>
    <div class="card-description">
        Visualizes the 'Big Picture' architecture of the complex.
    </div>
    <div class="card-features">
        <div class="card-features-title">Features:</div>
        <ul class="card-features-list">
            <li><strong>Outer Ring:</strong> Represents the full length of each selected protein sequence, proportionally sized by amino acid count</li>
            <li><strong>Inner Arcs:</strong> Curved lines connecting two genomic coordinates represent a crosslink between specific residue positions</li>
            <li><strong>Confidence Styles:</strong> Link opacity and thickness indicate confidence levels (Score or # PSMs), with thicker/opaque lines for higher confidence</li>
            <li><strong>Color Coding:</strong> Each protein has a unique color matching the Linear View palette for consistency</li>
            <li><strong>Domain Annotations:</strong> Functional domains are shown as colored blocks on the inner track, aligned with sequence positions</li>
        </ul>
    </div>
</div>
""", unsafe_allow_html=True)
