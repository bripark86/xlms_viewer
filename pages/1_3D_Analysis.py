import streamlit as st
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, Select
from io import StringIO
import streamlit.components.v1 as components
import json
import re
from difflib import get_close_matches
import os
import urllib.request

# Page configuration
st.set_page_config(page_title="3D Analysis", layout="wide")

# STRICT 6LTJ MAPPING (No paralogs allowed)
# Only maps proteins that physically exist in the 6LTJ PDB structure
BAF_MASTER_MAP = {
    # ATPase Module
    "SMARCA4": "I", "BRG1": "I", "SMCA4": "I",
    
    # Core Module
    "SMARCC2": "N", "BAF170": "N",
    "SMARCB1": "M", "SNF5": "M", "BAF47": "M",
    
    # Base Module
    "ARID1A": "L", "BAF250A": "L",
    
    # Actin Module
    "ACTB": "K", "Beta-actin": "K",
    "ACTL6A": "J", "BAF53A": "J",
    
    # BAF60
    "SMARCD1": "P", "BAF60A": "P",
    
    # BAF57
    "SMARCE1": "Q", "BAF57": "Q",
    
    # DPF2
    "DPF2": "R", "BAF45D": "R",
    
    # Histones (Nucleosome)
    "H3": ["A", "E"], "H3.3": ["A", "E"], "H33": ["A", "E"],
    "H4": ["B", "F"],
    "H2A": ["C", "G"],
    "H2B": ["D", "H"]
}

# File info mapping (same as 2D Analysis page)
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

# Initialize session state
if 'pdb_structure' not in st.session_state:
    st.session_state.pdb_structure = None
if 'pdb_filename' not in st.session_state:
    st.session_state.pdb_filename = None
if 'pdb_content' not in st.session_state:
    st.session_state.pdb_content = None
if 'pdb_chain_mapping' not in st.session_state:
    st.session_state.pdb_chain_mapping = {}
if 'csv_data' not in st.session_state:
    st.session_state.csv_data = None
if 'csv_data_processed' not in st.session_state:
    st.session_state.csv_data_processed = None
if 'protein_mapping' not in st.session_state:
    st.session_state.protein_mapping = {}
if 'auto_mapped' not in st.session_state:
    st.session_state.auto_mapped = {}
if 'links_data_for_table' not in st.session_state:
    st.session_state.links_data_for_table = []
if 'resolved_links_count' not in st.session_state:
    st.session_state.resolved_links_count = 0
if 'selected_dataset_3d' not in st.session_state:
    st.session_state.selected_dataset_3d = None

def download_pdb(pdb_id, filename):
    """Download PDB file from RCSB."""
    try:
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        urllib.request.urlretrieve(url, filename)
        return True
    except Exception as e:
        st.error(f"Error downloading PDB file: {str(e)}")
        return False

def load_default_pdb():
    """Load default 6LTJ.pdb, downloading if necessary."""
    default_file = "6LTJ.pdb"
    
    # Check if file exists locally
    if not os.path.exists(default_file):
        with st.spinner("Downloading default structure (6LTJ)..."):
            if download_pdb("6LTJ", default_file):
                st.success("Downloaded 6LTJ.pdb")
            else:
                return None, None
    
    # Load the file
    try:
        with open(default_file, 'r') as f:
            content = f.read()
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', StringIO(content))
        return structure, content
    except Exception as e:
        st.error(f"Error loading default PDB file: {str(e)}")
        return None, None

def parse_pdb_file(uploaded_file):
    """Parse PDB file and return structure object and raw content."""
    try:
        parser = PDBParser(QUIET=True)
        content = uploaded_file.read().decode('utf-8')
        structure = parser.get_structure('protein', StringIO(content))
        return structure, content
    except Exception as e:
        st.error(f"Error parsing PDB file: {str(e)}")
        return None, None

def get_pdb_chain_mapping(pdb_content, structure):
    """
    Extract chain mappings from PDB COMPND records.
    
    Parameters:
    -----------
    pdb_content : str
        Raw PDB file content
    structure : Bio.PDB.Structure
        Parsed structure object
    
    Returns:
    --------
    dict : Mapping of molecule names to list of chain IDs
           e.g., {'Histone H3': ['A', 'E'], 'SMARCB1': ['M']}
    """
    chain_mapping = {}
    
    # Parse COMPND records from header
    lines = pdb_content.split('\n')
    current_molecule = None
    current_chains = []
    
    for line in lines:
        if line.startswith('COMPND'):
            # COMPND format: COMPND   2 MOLECULE: Histone H3.1;
            # or COMPND   3 CHAIN: A, E;
            line_upper = line.upper()
            
            if 'MOLECULE:' in line_upper:
                # Extract molecule name
                mol_match = re.search(r'MOLECULE:\s*(.+?)(?:\s*;|$)', line, re.IGNORECASE)
                if mol_match:
                    # Save previous molecule if any
                    if current_molecule and current_chains:
                        chain_mapping[current_molecule] = list(set(current_chains))
                    
                    # Start new molecule
                    current_molecule = mol_match.group(1).strip().rstrip(';')
                    current_chains = []
            
            elif 'CHAIN:' in line_upper and current_molecule:
                # Extract chain IDs
                chain_match = re.search(r'CHAIN:\s*(.+?)(?:\s*;|$)', line, re.IGNORECASE)
                if chain_match:
                    chains_str = chain_match.group(1).strip().rstrip(';')
                    chains = [c.strip() for c in chains_str.split(',') if c.strip()]
                    current_chains.extend(chains)
    
    # Save last molecule
    if current_molecule and current_chains:
        chain_mapping[current_molecule] = list(set(current_chains))
    
    # Also add direct chain-to-name mapping from structure
    # (fallback if COMPND records are incomplete)
    model = structure[0]
    for chain in model:
        chain_id = chain.id
        # Try to infer from chain annotation or use chain ID as fallback
        if not any(chain_id in chains for chains in chain_mapping.values()):
            # If chain not found in COMPND, add it with a generic name
            chain_mapping[f'Chain {chain_id}'] = [chain_id]
    
    return chain_mapping

def get_ca_coordinates(structure, chain_id, resnum):
    """Get CA (Alpha Carbon) coordinates for a specific residue."""
    try:
        model = structure[0]  # Get first model
        if chain_id not in model:
            return None
        chain = model[chain_id]
        
        # Bio.PDB residue IDs are tuples: (hetero_flag, sequence_number, insertion_code)
        # Try to find residue by sequence number
        for residue in chain:
            # residue.id is a tuple: (hetero_flag, sequence_number, insertion_code)
            if residue.id[1] == resnum:  # Match sequence number
                if 'CA' in residue:
                    ca_atom = residue['CA']
                    coord = ca_atom.coord
                    return coord
        
        return None
    except Exception as e:
        return None

def calculate_shortest_distance(structure, res1_num, chains1_list, res2_num, chains2_list, protein1_name="", protein2_name=""):
    """
    Calculate shortest distance between two residues across multiple chain possibilities.
    
    Parameters:
    -----------
    structure : Bio.PDB.Structure
        Parsed PDB structure
    res1_num : int
        Residue number for first protein (will be cast to int)
    chains1_list : list
        List of chain IDs for first protein (e.g., ['A', 'E'])
    res2_num : int
        Residue number for second protein (will be cast to int)
    chains2_list : list
        List of chain IDs for second protein (e.g., ['M'])
    protein1_name : str
        Name of first protein (for debugging)
    protein2_name : str
        Name of second protein (for debugging)
    
    Returns:
    --------
    tuple : ((coord1, coord2, distance, chain1, chain2), []) on success, or (None, [missing_residues]) on failure
    """
    # Force integer casting
    res1_num = int(res1_num)
    res2_num = int(res2_num)
    
    min_distance = float('inf')
    best_coord1 = None
    best_coord2 = None
    best_chain1 = None
    best_chain2 = None
    missing_residues = []
    
    # Try all combinations of chains
    for chain1 in chains1_list:
        coord1 = get_ca_coordinates(structure, chain1.strip(), res1_num)
        if coord1 is None:
            missing_residues.append(f"Chain {chain1.strip()}, Res {res1_num}")
            continue
        
        for chain2 in chains2_list:
            coord2 = get_ca_coordinates(structure, chain2.strip(), res2_num)
            if coord2 is None:
                missing_residues.append(f"Chain {chain2.strip()}, Res {res2_num}")
                continue
            
            # Calculate Euclidean distance
            distance = np.linalg.norm(coord1 - coord2)
            
            if distance < min_distance:
                min_distance = distance
                best_coord1 = coord1
                best_coord2 = coord2
                best_chain1 = chain1.strip()
                best_chain2 = chain2.strip()
    
    # Return None with missing residues for debugging
    if best_coord1 is None or best_coord2 is None:
        return None, missing_residues
    
    # Convert numpy arrays to Python lists for JSON serialization
    coord1_list = best_coord1.tolist() if hasattr(best_coord1, 'tolist') else list(best_coord1)
    coord2_list = best_coord2.tolist() if hasattr(best_coord2, 'tolist') else list(best_coord2)
    
    return (coord1_list, coord2_list, min_distance, best_chain1, best_chain2), []

def clean_protein_name(name):
    """Clean protein name by removing prefix (e.g., 'sp|Q96019|ACTL6A' -> 'ACTL6A')."""
    if pd.isna(name) or name == "":
        return name
    if '|' in str(name):
        parts = str(name).split('|')
        return parts[-1]
    return str(name)

def parse_csv_file(uploaded_file):
    """Parse CSV file and standardize column names."""
    try:
        df = pd.read_csv(uploaded_file)
        
        # Detect format and standardize
        if 'PROTEIN.1' in df.columns and 'PROTEIN.2' in df.columns:
            # Proxl format
            df = df.copy()
            df['Protein1'] = df['PROTEIN.1']
            df['LinkPos1'] = pd.to_numeric(df['POSITION'], errors='coerce')
            df['Protein2'] = df['PROTEIN.2']
            df['LinkPos2'] = pd.to_numeric(df['POSITION.1'], errors='coerce')
        elif 'Protein1' not in df.columns:
            # Try case-insensitive matching
            cols_lower = {col.lower(): col for col in df.columns}
            if 'protein1' in cols_lower and 'protein2' in cols_lower:
                df['Protein1'] = df[cols_lower['protein1']]
                df['Protein2'] = df[cols_lower['protein2']]
                if 'linkpos1' in cols_lower:
                    df['LinkPos1'] = pd.to_numeric(df[cols_lower['linkpos1']], errors='coerce')
                if 'linkpos2' in cols_lower:
                    df['LinkPos2'] = pd.to_numeric(df[cols_lower['linkpos2']], errors='coerce')
        
        # Clean protein names
        if 'Protein1' in df.columns:
            df['Protein1_clean'] = df['Protein1'].apply(clean_protein_name)
            df['Protein2_clean'] = df['Protein2'].apply(clean_protein_name)
        
        return df
    except Exception as e:
        st.error(f"Error parsing CSV file: {str(e)}")
        return None

def load_dataset_3d(dataset_key):
    """Load dataset from data/ folder for 3D analysis."""
    if dataset_key not in FILE_INFO_LIST:
        return None
    
    info = FILE_INFO_LIST[dataset_key]
    stem = info['stem']
    dataset_type = info['type']
    
    try:
        if dataset_type == "standard":
            links_path = f"data/{stem}_Links.csv"
            if not os.path.exists(links_path):
                return None
            df = pd.read_csv(links_path)
            
        elif dataset_type == "standard_shortname":
            links_path = f"data/{stem}.csv"
            if not os.path.exists(links_path):
                return None
            df = pd.read_csv(links_path)
            
            # Map short names using protein_lengths if available
            try:
                protein_lengths = pd.read_csv("data/protein_lengths.csv")
                name_to_rawid = protein_lengths[['short_name', 'raw_id']].drop_duplicates()
                
                df = df.merge(
                    name_to_rawid, left_on='Protein1', right_on='short_name', how='left'
                ).rename(columns={'raw_id': 'Protein1_raw'}).merge(
                    name_to_rawid, left_on='Protein2', right_on='short_name', how='left'
                ).rename(columns={'raw_id': 'Protein2_raw'})
                
                df = df.dropna(subset=['Protein1_raw', 'Protein2_raw'])
                df = df.drop(columns=['Protein1', 'Protein2'])
                df = df.rename(columns={'Protein1_raw': 'Protein1', 'Protein2_raw': 'Protein2'})
            except Exception:
                # If protein_lengths.csv doesn't exist, use as-is
                pass
        
        else:
            return None
        
        # Clean protein names
        if 'Protein1' in df.columns:
            df['Protein1_clean'] = df['Protein1'].apply(clean_protein_name)
            df['Protein2_clean'] = df['Protein2'].apply(clean_protein_name)
        
        return df
    
    except Exception as e:
        st.error(f"Error loading dataset {dataset_key}: {str(e)}")
        return None

def auto_apply_baf_master_map(df):
    """Automatically apply BAF_MASTER_MAP to protein mapping."""
    if df is None or df.empty:
        return {}
    
    mapping = {}
    unique_proteins = set()
    if 'Protein1_clean' in df.columns:
        unique_proteins.update(df['Protein1_clean'].dropna().unique())
        unique_proteins.update(df['Protein2_clean'].dropna().unique())
    
    for protein in unique_proteins:
        protein_upper = protein.upper().strip()
        
        # Check BAF_MASTER_MAP
        for alias, chains in BAF_MASTER_MAP.items():
            if protein_upper == alias.upper():
                # Handle both string and list formats
                if isinstance(chains, list):
                    mapping[protein] = ', '.join(sorted(chains))
                else:
                    mapping[protein] = str(chains)
                break
    
    return mapping

def filter_mapped_proteins(df, protein_mapping):
    """
    Filter out crosslinks where proteins don't have valid chain mappings.
    
    Returns filtered dataframe and count of filtered rows.
    """
    if df is None or df.empty:
        return df, 0
    
    original_count = len(df)
    
    # Filter rows where both proteins have valid mappings
    mask = df.apply(lambda row: (
        pd.notna(row.get('Protein1_clean')) and
        pd.notna(row.get('Protein2_clean')) and
        str(row['Protein1_clean']) in protein_mapping and
        str(row['Protein2_clean']) in protein_mapping and
        protein_mapping[str(row['Protein1_clean'])].strip() != '' and
        protein_mapping[str(row['Protein2_clean'])].strip() != ''
    ), axis=1)
    
    filtered_df = df[mask].copy()
    filtered_count = original_count - len(filtered_df)
    
    return filtered_df, filtered_count

def generate_ngl_html(pdb_content, links_data, focus_target=None):
    """
    Generate HTML string with embedded NGL.js viewer.
    
    Parameters:
    -----------
    pdb_content : str
        PDB file content as string
    links_data : list
        List of dictionaries with keys: 'start' (list[x,y,z]), 'end' (list[x,y,z]), 
        'color' (str hex), 'radius' (float)
    focus_target : list, optional
        List of [x, y, z] coordinates to zoom to. If None, auto-view is used.
    
    Returns:
    --------
    str : Complete HTML string with embedded NGL.js viewer
    """
    # Escape the PDB content for JavaScript
    pdb_json = json.dumps(pdb_content)
    
    # Convert links_data to JSON
    links_json = json.dumps(links_data)
    
    # Convert focus_target to JSON (if provided)
    focus_target_json = json.dumps(focus_target) if focus_target else "null"
    
    html_template = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <script src="https://unpkg.com/ngl@2.0.0-dev.37/dist/ngl.js"></script>
        <style>
            body {{
                margin: 0;
                padding: 0;
                overflow: hidden;
                background-color: #000000;
            }}
            #viewport {{
                width: 100%;
                height: 600px;
            }}
        </style>
    </head>
    <body>
        <div id="viewport"></div>
        <script>
            // PDB content
            const pdbContent = {pdb_json};
            const linksData = {links_json};
            const focusTarget = {focus_target_json};
            
            // Initialize NGL Stage
            const stage = new NGL.Stage("viewport", {{
                backgroundColor: "black"
            }});
            
            // Load PDB structure from blob
            const blob = new Blob([pdbContent], {{ type: "text/plain" }});
            const file = new File([blob], "structure.pdb", {{ type: "text/plain" }});
            
            stage.loadFile(file, {{ ext: "pdb" }}).then(function (component) {{
                // Add cartoon representation with chainindex color scheme
                component.addRepresentation("cartoon", {{
                    colorScheme: "chainindex"
                }});
                
                // Add crosslinks as cylinders
                if (linksData && linksData.length > 0) {{
                    var shape = new NGL.Shape("crosslinks");
                    
                    linksData.forEach(function(link) {{
                        // Create NGL Vectors explicitly to ensure valid coords
                        var s = new NGL.Vector3(link.start[0], link.start[1], link.start[2]);
                        var e = new NGL.Vector3(link.end[0], link.end[1], link.end[2]);
                        
                        // Convert hex color to RGB array (default to yellow [1, 1, 0])
                        var colorArray = [1, 1, 0]; // Yellow default
                        if (link.color) {{
                            // Parse hex color (e.g., "#FFFF00" -> [1, 1, 0])
                            var hex = link.color.replace('#', '');
                            if (hex.length === 6) {{
                                colorArray = [
                                    parseInt(hex.substr(0, 2), 16) / 255,
                                    parseInt(hex.substr(2, 2), 16) / 255,
                                    parseInt(hex.substr(4, 2), 16) / 255
                                ];
                            }}
                        }}
                        
                        shape.addCylinder(s, e, colorArray, link.radius);
                    }});
                    
                    var shapeComp = stage.addComponentFromObject(shape);
                    shapeComp.addRepresentation("buffer");
                }}
                
                // Auto-zoom to fit
                if (focusTarget) {{
                    // Zoom to specific link
                    var focusVec = new NGL.Vector3(focusTarget[0], focusTarget[1], focusTarget[2]);
                    stage.animationControls.zoom(focusVec, 2000);
                }} else {{
                    // Default auto view
                    stage.autoView(2000);
                }}
            }});
        </script>
    </body>
    </html>
    """
    return html_template

# Main UI
st.title("🔬 3D Analysis")
st.markdown("Visualize cross-linking mass spectrometry data on 3D protein structures")

# Sidebar
with st.sidebar:
    st.header("📁 Input Files")
    
    # Dataset selector
    dataset_options = {
        "Upload Custom CSV": None,
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
        "Select Dataset",
        options=list(dataset_options.keys()),
        key="dataset_selector_3d"
    )
    
    selected_dataset_key = dataset_options[selected_dataset_display]
    
    # Load dataset if selected
    if selected_dataset_key is not None:
        if st.session_state.selected_dataset_3d != selected_dataset_key:
            with st.spinner(f"Loading dataset: {selected_dataset_display}..."):
                df = load_dataset_3d(selected_dataset_key)
                if df is not None:
                    st.session_state.csv_data = df
                    st.session_state.selected_dataset_3d = selected_dataset_key
                    st.success(f"Loaded: {selected_dataset_display} ({len(df)} links)")
                    
                    # Auto-apply BAF_MASTER_MAP
                    auto_mapping = auto_apply_baf_master_map(df)
                    if auto_mapping:
                        st.session_state.protein_mapping = auto_mapping
                        st.success(f"✅ Auto-mapped {len(auto_mapping)} protein(s) using BAF_MASTER_MAP")
                    
                    # Process and filter immediately if PDB is loaded
                    if st.session_state.pdb_structure is not None and auto_mapping:
                        filtered_df, filtered_count = filter_mapped_proteins(df, auto_mapping)
                        st.session_state.csv_data_processed = filtered_df
                        st.session_state.resolved_links_count = 0
                        if filtered_count > 0:
                            st.warning(f"⚠️ {filtered_count} link(s) hidden (missing chain mapping)")
                else:
                    st.error(f"Could not load dataset: {selected_dataset_display}")
    
    # CSV file upload (only show if "Upload Custom CSV" is selected)
    csv_file = None
    if selected_dataset_key is None:
        csv_file = st.file_uploader("Upload Crosslinks CSV", type=['csv'], key="csv_upload")
        if csv_file is not None:
            with st.spinner("Parsing CSV file..."):
                df = parse_csv_file(csv_file)
                if df is not None:
                    st.session_state.csv_data = df
                    st.session_state.selected_dataset_3d = None  # Reset dataset selection
                    st.session_state.csv_data_processed = None  # Reset processed data
                    st.session_state.protein_mapping = {}  # Reset mappings
                    st.success(f"Loaded: {csv_file.name} ({len(df)} links)")
    
    # PDB file upload
    pdb_file = st.file_uploader("Upload PDB File", type=['pdb', 'ent'], key="pdb_upload")
    
    # Load default PDB if no file uploaded
    if pdb_file is None:
        if st.session_state.pdb_structure is None or st.session_state.pdb_filename != "6LTJ.pdb":
            with st.spinner("Loading default structure (6LTJ)..."):
                structure, content = load_default_pdb()
                if structure is not None and content is not None:
                    st.session_state.pdb_structure = structure
                    st.session_state.pdb_content = content
                    st.session_state.pdb_filename = "6LTJ.pdb"
                    
                    # Extract chain mappings from PDB header
                    chain_mapping = get_pdb_chain_mapping(content, structure)
                    st.session_state.pdb_chain_mapping = chain_mapping
                    
                    st.info("📦 Using default structure: 6LTJ.pdb")
                    if chain_mapping:
                        st.caption(f"Found {len(chain_mapping)} molecule(s) in PDB header")
                    
                    # Reset mapping when new PDB is loaded
                    st.session_state.protein_mapping = {}
                    st.session_state.auto_mapped = {}
                    st.session_state.csv_data_processed = None
                    st.session_state.resolved_links_count = 0
    else:
        if st.session_state.pdb_filename != pdb_file.name:
            with st.spinner("Parsing PDB file and extracting metadata..."):
                structure, content = parse_pdb_file(pdb_file)
                if structure is not None and content is not None:
                    st.session_state.pdb_structure = structure
                    st.session_state.pdb_content = content
                    st.session_state.pdb_filename = pdb_file.name
                    
                    # Extract chain mappings from PDB header
                    chain_mapping = get_pdb_chain_mapping(content, structure)
                    st.session_state.pdb_chain_mapping = chain_mapping
                    
                    st.success(f"Loaded: {pdb_file.name}")
                    if chain_mapping:
                        st.caption(f"Found {len(chain_mapping)} molecule(s) in PDB header")
                    
                    # Reset mapping when new PDB is loaded
                    st.session_state.protein_mapping = {}
                    st.session_state.auto_mapped = {}
                    st.session_state.csv_data_processed = None
                    st.session_state.resolved_links_count = 0
    
    # Process CSV data if available (only if not already processed from dataset loading)
    if st.session_state.csv_data is not None and st.session_state.csv_data_processed is None:
        df = st.session_state.csv_data
        
        # Auto-map proteins if PDB is loaded and mapping not already set
        if st.session_state.pdb_structure is not None and st.session_state.pdb_chain_mapping:
            if not st.session_state.protein_mapping:
                unique_proteins = set()
                if 'Protein1_clean' in df.columns:
                    unique_proteins.update(df['Protein1_clean'].dropna().unique())
                    unique_proteins.update(df['Protein2_clean'].dropna().unique())
                unique_proteins = sorted([p for p in unique_proteins if p])
                
                if unique_proteins:
                    # Try BAF_MASTER_MAP first
                    auto_mapping = auto_apply_baf_master_map(df)
                    if not auto_mapping:
                        # Fallback: try to match with PDB chain mapping
                        pdb_molecules = list(st.session_state.pdb_chain_mapping.keys())
                        for csv_protein in unique_proteins:
                            csv_protein_upper = csv_protein.upper().strip()
                            # Try exact match in master map
                            matched_chains = None
                            for alias, chains in BAF_MASTER_MAP.items():
                                if csv_protein_upper == alias.upper():
                                    if isinstance(chains, list):
                                        matched_chains = chains
                                    else:
                                        matched_chains = [chains]
                                    auto_mapping[csv_protein] = ', '.join(sorted(matched_chains))
                                    break
                    
                    st.session_state.protein_mapping = auto_mapping
                    if auto_mapping:
                        st.success(f"✅ Auto-mapped {len(auto_mapping)} protein(s)")
        
        # Process and filter the data
        if st.session_state.protein_mapping:
            filtered_df, filtered_count = filter_mapped_proteins(df, st.session_state.protein_mapping)
            st.session_state.csv_data_processed = filtered_df
            st.session_state.resolved_links_count = 0  # Reset until recalculated
            
            if filtered_count > 0:
                st.warning(f"⚠️ {filtered_count} link(s) hidden (missing chain mapping)")
    
    st.divider()
    
    st.header("🎨 Visualization Settings")
    link_color = st.color_picker("Link Color", value="#FFFF00", key="link_color")
    link_radius = 0.20  # Fixed link radius
    
    st.divider()
    
    # Protein-Chain Mapping Section
    st.header("🔗 Protein-Chain Mapping")
    
    if st.session_state.csv_data is None:
        st.info("Select a dataset or upload a CSV file to see proteins")
    else:
        df = st.session_state.csv_data
        
        if 'Protein1_clean' not in df.columns:
            st.error("CSV format not recognized")
        else:
            # Get unique proteins from CSV
            unique_proteins_csv = set()
            unique_proteins_csv.update(df['Protein1_clean'].dropna().unique())
            unique_proteins_csv.update(df['Protein2_clean'].dropna().unique())
            unique_proteins_csv = sorted([p for p in unique_proteins_csv if p])
            
            # Filter to only show proteins that exist in the PDB structure
            # Get all chains that exist in the PDB
            pdb_chains = set()
            if st.session_state.pdb_chain_mapping:
                for chains_list in st.session_state.pdb_chain_mapping.values():
                    pdb_chains.update(chains_list)
            
            # Also get chains from the structure object directly
            if st.session_state.pdb_structure:
                for model in st.session_state.pdb_structure:
                    for chain in model:
                        pdb_chains.add(chain.id)
            
            # Filter proteins: only show those that have mappings to existing PDB chains
            unique_proteins = []
            for protein in unique_proteins_csv:
                # Check if protein has a mapping
                has_valid_mapping = False
                
                # Check session state mapping first
                if protein in st.session_state.protein_mapping:
                    mapping = st.session_state.protein_mapping[protein].strip()
                    if mapping:
                        mapped_chains = [c.strip() for c in mapping.split(',') if c.strip()]
                        if any(chain in pdb_chains for chain in mapped_chains):
                            has_valid_mapping = True
                
                # Check BAF_MASTER_MAP if no session state mapping
                if not has_valid_mapping and protein in BAF_MASTER_MAP:
                    map_value = BAF_MASTER_MAP[protein]
                    if isinstance(map_value, list):
                        if any(chain in pdb_chains for chain in map_value):
                            has_valid_mapping = True
                    elif map_value in pdb_chains:
                        has_valid_mapping = True
                
                if has_valid_mapping:
                    unique_proteins.append(protein)
            
            if not unique_proteins:
                st.info("No proteins from the CSV match chains in the loaded PDB structure.")
            else:
                st.markdown("**Map each protein to chain ID(s). Use commas for multiple chains (e.g., `A, E`).**")
                if len(unique_proteins) < len(unique_proteins_csv):
                    st.markdown(f"*Showing {len(unique_proteins)} of {len(unique_proteins_csv)} proteins that exist in the PDB structure.*")
                
                # Show PDB metadata if available
                if st.session_state.pdb_chain_mapping:
                    with st.expander("📋 PDB Metadata (COMPND Records)", expanded=False):
                        # Create a filtered and renamed mapping for display
                        display_mapping = {}
                        for mol_name, chains in st.session_state.pdb_chain_mapping.items():
                            # Skip Chain X - don't show it
                            chains_list = chains if isinstance(chains, list) else [chains]
                            if 'X' in chains_list:
                                continue
                            if mol_name == 'Chain X':
                                continue
                            
                            # Rename specific chains with full protein names and aliases
                            mol_name_upper = mol_name.upper()
                            chains_list = chains if isinstance(chains, list) else [chains]
                            
                            if 'M' in chains_list and len(chains_list) == 1 and ('Chain M' in mol_name or mol_name == 'M'):
                                display_mapping['SWI/SNF-related matrix-associated actin-dependent regulator of chromatin subfamily B member(SMARCB1)'] = chains
                            elif 'P' in chains_list and len(chains_list) == 1 and ('Chain P' in mol_name or mol_name == 'P'):
                                display_mapping['SWI/SNF-related matrix-associated actin-dependent regulator of chromatin subfamily D member 1(SMARCD1)'] = chains
                            elif 'TRANSCRIPTION ACTIVATOR BRG1' in mol_name_upper and 'I' in chains_list:
                                display_mapping['TRANSCRIPTION ACTIVATOR BRG1(SMARCA4)'] = chains
                            elif 'ZINC FINGER PROTEIN UBI-D4' in mol_name_upper and 'R' in chains_list:
                                display_mapping['ZINC FINGER PROTEIN UBI-D4(DPF2)'] = chains
                            elif 'SWI/SNF-RELATED MATRIX-ASSOCIATED ACTIN-DEPENDENT REGULATOR' in mol_name_upper and 'Q' in chains_list and 'SMARCE1' not in mol_name_upper:
                                display_mapping['SWI/SNF-RELATED MATRIX-ASSOCIATED ACTIN-DEPENDENT REGULATOR(SMARCE1)'] = chains
                            elif 'ACTIN, CYTOPLASMIC 1' in mol_name_upper and 'K' in chains_list:
                                display_mapping['ACTIN, CYTOPLASMIC 1(ACTB)'] = chains
                            elif 'AT-RICH INTERACTIVE DOMAIN-CONTAINING PROTEIN 1A' in mol_name_upper and 'L' in chains_list:
                                display_mapping['AT-RICH INTERACTIVE DOMAIN-CONTAINING PROTEIN 1A(ARID1A)'] = chains
                            elif 'ACTIN-LIKE PROTEIN 6A' in mol_name_upper and 'J' in chains_list:
                                display_mapping['ACTIN-LIKE PROTEIN 6A(ACTL6A)'] = chains
                            else:
                                display_mapping[mol_name] = chains
                        
                        for mol_name, chains in sorted(display_mapping.items()):
                            chains_list = chains if isinstance(chains, list) else [chains]
                            st.text(f"{mol_name}: {', '.join(sorted(chains_list))}")
                
                st.markdown("")
                
                # Create mapping inputs (without auto-mapping indicators)
                for protein in unique_proteins:
                    current_value = st.session_state.protein_mapping.get(protein, "")
                    
                    label = f"**{protein}** → Chain(s)"
                    help_text = f"Enter chain ID(s) for {protein}, separated by commas if multiple"
                    
                    new_value = st.text_input(
                        label,
                        value=current_value,
                        key=f"mapping_{protein}",
                        help=help_text
                    )
                    st.session_state.protein_mapping[protein] = new_value

# Main content area
col1, col2 = st.columns([2, 1])

with col1:
    st.header("3D Structure View")
    
    if st.session_state.pdb_structure is None:
        st.info("👈 Please upload a PDB file to begin")
    else:
        # Get PDB content for rendering
        pdb_content = st.session_state.pdb_content
        if pdb_content is not None:
            # Build links data for NGL.js
            links_data = []
            links_added = 0
            links_failed = 0
            all_missing_residues = []
            
            if st.session_state.csv_data_processed is not None and st.session_state.protein_mapping:
                df_processed = st.session_state.csv_data_processed.copy()
                structure = st.session_state.pdb_structure
                
                # Render ALL valid crosslinks (no selection filter)
                df_to_render = df_processed.copy()
                
                # Track statistics
                mapped_count = 0
                resolved_count = 0
                
                with st.spinner("Calculating crosslinks..."):
                    links_data = []
                    link_details = []  # Store details for Link Inspector table
                    all_missing_residues = []
                    links_failed = 0
                    links_added = 0
                    mapped_count = 0
                    resolved_count = 0
                    
                    for idx, row in df_to_render.iterrows():
                        if pd.isna(row.get('Protein1_clean')) or pd.isna(row.get('Protein2_clean')):
                            continue
                        if pd.isna(row.get('LinkPos1')) or pd.isna(row.get('LinkPos2')):
                            continue
                        
                        protein1 = str(row['Protein1_clean'])
                        protein2 = str(row['Protein2_clean'])
                        
                        # Force integer casting (CSV may have floats like 100.0)
                        try:
                            res1 = int(float(row['LinkPos1']))
                            res2 = int(float(row['LinkPos2']))
                        except (ValueError, TypeError):
                            links_failed += 1
                            continue
                        
                        # Get chain mappings
                        chains1_str = st.session_state.protein_mapping.get(protein1, "")
                        chains2_str = st.session_state.protein_mapping.get(protein2, "")
                        
                        if not chains1_str or not chains2_str:
                            links_failed += 1
                            continue
                        
                        # Count as mapped (has valid protein mapping)
                        mapped_count += 1
                        
                        chains1 = [c.strip() for c in chains1_str.split(',') if c.strip()]
                        chains2 = [c.strip() for c in chains2_str.split(',') if c.strip()]
                        
                        if not chains1 or not chains2:
                            links_failed += 1
                            continue
                        
                        # Calculate shortest distance (with debugging)
                        result, missing = calculate_shortest_distance(structure, res1, chains1, res2, chains2, protein1, protein2)
                        
                        if result is not None:
                            # Success: result is (coord1, coord2, distance, chain1, chain2)
                            coord1, coord2, distance, chain1, chain2 = result
                            resolved_count += 1
                            
                            # Add to links_data for NGL.js
                            # coord1 and coord2 are already Python lists (converted from numpy arrays)
                            links_data.append({
                                'start': coord1,  # Already a list
                                'end': coord2,    # Already a list
                                'color': link_color,
                                'radius': link_radius
                            })
                            
                            # Store details for Link Inspector table
                            link_details.append({
                                'Protein1': protein1,
                                'Res1': res1,
                                'Chain1': chain1,
                                'Protein2': protein2,
                                'Res2': res2,
                                'Chain2': chain2,
                                'Distance': distance,
                                'Score': row.get('Score', 'N/A') if 'Score' in row else 'N/A',
                                'midpoint': [(coord1[0] + coord2[0]) / 2, (coord1[1] + coord2[1]) / 2, (coord1[2] + coord2[2]) / 2],
                                'index': len(link_details)  # Index in links_data array
                            })
                            links_added += 1
                        else:
                            # Failure: missing contains list of missing residues
                            if missing:
                                all_missing_residues.extend(missing)
                            links_failed += 1
                    
                    # Missing residue warnings removed - no longer displayed
                    
                    # Store link details for Link Inspector
                    st.session_state.links_data_for_table = link_details
                    
                    # Store the count of successfully resolved links for use in Crosslink table
                    st.session_state.resolved_links_count = links_added
                    
                    # Calculate and display detailed status
                    missing_mapping_count = len(df_to_render) - mapped_count
                    missing_residue_count = mapped_count - resolved_count
                    
                    # Status message removed as requested - no longer displayed above 3D view
            
            # Generate HTML with NGL.js viewer (no zoom functionality)
            html_content = generate_ngl_html(pdb_content, links_data, focus_target=None)
            
            
            # Render the HTML component
            components.html(html_content, height=600, width=800)
            
            # Link Inspector table
            if st.session_state.links_data_for_table and len(st.session_state.links_data_for_table) > 0:
                st.divider()
                st.header("🔍 Link Inspector")
                st.caption("Crosslink details for currently visualized links")
                
                # Create dataframe for display
                inspector_data = []
                for i, link in enumerate(st.session_state.links_data_for_table):
                    inspector_data.append({
                        'Cylinder #': i,  # 0-indexed to match 3D viewer display
                        'Protein1': link['Protein1'],
                        'Res1': link['Res1'],
                        'Chain1': link['Chain1'],
                        'Protein2': link['Protein2'],
                        'Res2': link['Res2'],
                        'Chain2': link['Chain2'],
                        'Distance (Å)': f"{link['Distance']:.1f}",
                        'Score': link['Score']
                    })
                
                inspector_df = pd.DataFrame(inspector_data)
                
                # Display as static dataframe (no selection/zoom functionality)
                st.dataframe(
                    inspector_df,
                    width='stretch',
                    hide_index=True
                )

with col2:
    st.header("📊 Crosslink Table")
    
    # Interactive crosslink table with "Show" checkbox
    if st.session_state.csv_data_processed is not None and st.session_state.protein_mapping:
        df_processed = st.session_state.csv_data_processed.copy()
        
        # Filter table to ONLY show links that have valid mappings (can be drawn)
        if 'Protein1_clean' in df_processed.columns and 'Protein2_clean' in df_processed.columns:
            valid_mask = df_processed.apply(lambda row: (
                pd.notna(row.get('Protein1_clean')) and
                pd.notna(row.get('Protein2_clean')) and
                str(row['Protein1_clean']) in st.session_state.protein_mapping and
                str(row['Protein2_clean']) in st.session_state.protein_mapping and
                st.session_state.protein_mapping[str(row['Protein1_clean'])].strip() != '' and
                st.session_state.protein_mapping[str(row['Protein2_clean'])].strip() != ''
            ), axis=1)
            df_processed = df_processed[valid_mask].copy()
        
        if len(df_processed) == 0:
            st.info("No crosslinks have valid protein-chain mappings. Check your mappings in the sidebar.")
        else:
            # Prepare display columns (no 'Show' column needed - all links are automatically visualized)
            display_cols = []
            if 'Protein1_clean' in df_processed.columns:
                display_cols.append('Protein1_clean')
            if 'LinkPos1' in df_processed.columns:
                display_cols.append('LinkPos1')
            if 'Protein2_clean' in df_processed.columns:
                display_cols.append('Protein2_clean')
            if 'LinkPos2' in df_processed.columns:
                display_cols.append('LinkPos2')
            if 'Score' in df_processed.columns:
                display_cols.append('Score')
            
            # Select only available columns
            display_cols = [col for col in display_cols if col in df_processed.columns]
            df_display = df_processed[display_cols].copy()
            
            # Rename for display
            rename_dict = {
                'Protein1_clean': 'Protein1',
                'LinkPos1': 'Pos1',
                'Protein2_clean': 'Protein2',
                'LinkPos2': 'Pos2'
            }
            df_display = df_display.rename(columns=rename_dict)
            
            # Display as simple dataframe (read-only, all links are automatically visualized)
            st.dataframe(
                df_display,
                width='stretch',
                hide_index=True
            )
            
            # Show statistics
            total_mapped = len(df_processed)
            total_loaded = len(st.session_state.csv_data) if st.session_state.csv_data is not None else 0
            total_unmapped = total_loaded - total_mapped
            
            # Use the actual count of successfully resolved links (from 3D visualization)
            resolved_count = st.session_state.resolved_links_count if st.session_state.resolved_links_count > 0 else 0
            status_msg = f"Visualizing {resolved_count} links"
            
            # Show breakdown of why some links aren't visualized
            if total_unmapped > 0:
                status_msg += f" | {total_unmapped} hidden (no mapping)"
            if total_mapped > resolved_count:
                missing_residue_count = total_mapped - resolved_count
                status_msg += f" | {missing_residue_count} hidden (residue not found in structure)"
            st.caption(status_msg)
    else:
        st.info("Complete the protein-chain mapping (in sidebar) to see crosslink selection table")

# Footer
st.divider()
st.markdown("**Usage Tips:**")
st.markdown("""
- The app loads **6LTJ.pdb** by default (downloads if needed)
- Select a dataset from the dropdown or upload your own crosslinks CSV file
- When using embedded datasets, protein mappings are automatically applied using BAF_MASTER_MAP
- Map each protein to its chain ID(s) in the sidebar if needed
- Use comma-separated values for proteins with multiple copies (e.g., H3 → `A, E`)
- **Hover over the structure** to see residue information
- All valid crosslinks are automatically visualized
- The app automatically selects the closest chain pair for each crosslink
- Adjust link color and radius in the sidebar
""")

