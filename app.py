import streamlit as st
import pandas as pd
from io import StringIO
from farsa_parser import FastaParser
from gene_classes import MitochondrialDNA, BioPythonAligner
import plotly.express as px
import plotly.graph_objects as go

st.set_page_config(page_title="DNA Analysis", layout="wide")

def main():
    st.title("DNA Analysis Tool")
    
    # Sidebar for navigation
    page = st.sidebar.selectbox(
        "Select Analysis",
        ["Upload & Overview", "Sequence Analysis", "Motif Analysis", "Sequence Alignment"]
    )
    
    # File uploader that persists across pages
    if 'data' not in st.session_state:
        st.session_state.data = None
        st.session_state.mtdna_objects = None
    
    uploaded_file = st.sidebar.file_uploader("Upload FASTA file", type=['txt', 'fasta'])
    if uploaded_file is not None:
        try:
            stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
            with open("temp_fasta.txt", "w") as f:
                f.write(stringio.read())
            
            # Parse the file and create MitochondrialDNA objects
            st.session_state.data = FastaParser.parse_file("temp_fasta.txt")
            st.session_state.mtdna_objects = FastaParser.create_mitochondrial_objects(st.session_state.data)
            
        except Exception as e:
            st.error(f"Error parsing file: {str(e)}")
    
    # Display different pages based on selection
    if page == "Upload & Overview":
        show_overview()
    elif page == "Sequence Analysis":
        show_sequence_analysis()
    elif page == "Motif Analysis":
        show_motif_analysis()
    elif page == "Sequence Alignment":
        show_sequence_alignment()

def show_overview():
    st.header("Dataset Overview")
    
    if st.session_state.data is not None:
        # Display basic statistics
        st.subheader("Dataset Summary")
        col1, col2 = st.columns(2)
        
        with col1:
            st.metric("Number of Sequences", len(st.session_state.data))
            
        with col2:
            avg_length = sum(len(obj.sequence) for obj in st.session_state.mtdna_objects) / len(st.session_state.mtdna_objects)
            st.metric("Average Sequence Length", f"{avg_length:,.0f}")
        
        # Display sequence lengths distribution
        lengths = [obj.sequence_length() for obj in st.session_state.mtdna_objects]
        fig = px.histogram(
            x=lengths,
            title="Distribution of Sequence Lengths",
            labels={'x': 'Sequence Length', 'y': 'Count'}
        )
        st.plotly_chart(fig)
        
        # Display data table
        st.subheader("Sequences in Dataset")
        df_display = st.session_state.data.copy()
        df_display['Sequence Length'] = lengths
        st.dataframe(df_display)
    else:
        st.info("Please upload a FASTA file to begin analysis")

def show_sequence_analysis():
    st.header("Sequence Analysis")
    
    if st.session_state.mtdna_objects:
        # Select sequence for analysis
        selected_seq = st.selectbox(
            "Select Sequence to Analyze",
            [obj.sequence_id for obj in st.session_state.mtdna_objects]
        )
        
        mtdna = next(obj for obj in st.session_state.mtdna_objects if obj.sequence_id == selected_seq)
        
        # Display sequence stats
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Sequence Length", mtdna.sequence_length())
        
        with col2:
            st.metric("GC Content", mtdna.GC_content_percentage())
        
        # Base composition analysis
        st.subheader("Base Composition")
        base_comp = mtdna.base_composition()
        base_comp_pct = mtdna.base_composition_percentage()
        
        # Create two column layout for absolute and percentage charts
        col1, col2 = st.columns(2)
        
        with col1:
            fig_abs = go.Figure(data=[
                go.Bar(
                    x=list(base_comp.keys()),
                    y=list(base_comp.values()),
                    marker_color=['#FF9999', '#66B2FF', '#99FF99', '#FFCC99']
                )
            ])
            fig_abs.update_layout(
                title="Absolute Base Composition",
                xaxis_title="Base",
                yaxis_title="Count",
                bargap=0.3,
                height=400,
                yaxis=dict(
                    range=[0, max(base_comp.values()) * 1.1]  # Add 10% padding
                )
            )
            st.plotly_chart(fig_abs, use_container_width=True)
            
        with col2:
            fig_pct = go.Figure(data=[
                go.Bar(
                    x=list(base_comp_pct.keys()),
                    y=list(base_comp_pct.values()),
                    marker_color=['#FF9999', '#66B2FF', '#99FF99', '#FFCC99']
                )
            ])
            fig_pct.update_layout(
                title="Relative Base Composition (%)",
                xaxis_title="Base",
                yaxis_title="Percentage",
                bargap=0.3,
                height=400,
                yaxis=dict(
                    range=[0, 100]  # Fix range for percentages
                )
            )
            st.plotly_chart(fig_pct, use_container_width=True)
        
        # Sequence viewer
        st.subheader("Sequence Viewer")
        sequence_window = st.slider(
            "Select sequence window",
            0,
            mtdna.sequence_length(),
            (0, min(1000, mtdna.sequence_length())),
            step=100
        )
        st.code(mtdna.extract_subsequences(sequence_window[0], sequence_window[1]))
        
    else:
        st.info("Please upload a FASTA file to begin analysis")

def show_motif_analysis():
    st.header("Motif Analysis")
    
    if st.session_state.mtdna_objects:
        col1, col2 = st.columns([1, 2])
        
        with col1:
            # Motif input
            motif = st.text_input(
                "Enter motif sequence to search",
                value="ATCG",
                help="Enter a DNA sequence motif to search for"
            ).upper()
            
            # Select sequences for analysis
            selected_seqs = st.multiselect(
                "Select Sequences to Analyze",
                [obj.sequence_id for obj in st.session_state.mtdna_objects],
                default=[obj.sequence_id for obj in st.session_state.mtdna_objects][:2]
            )
        
        if motif and selected_seqs:
            # Analyze motif occurrences
            results = []
            for seq_id in selected_seqs:
                mtdna = next(obj for obj in st.session_state.mtdna_objects if obj.sequence_id == seq_id)
                positions = mtdna.find_motif(motif)
                results.append({
                    'sequence_id': seq_id,
                    'occurrences': len(positions),
                    'positions': positions
                })
            
            # Display results
            with col2:
                st.subheader("Motif Analysis Results")
                
                # Bar chart of occurrences
                fig = px.bar(
                    results,
                    x='sequence_id',
                    y='occurrences',
                    title=f"Occurrences of motif '{motif}'",
                    labels={'sequence_id': 'Sequence', 'occurrences': 'Number of Occurrences'}
                )
                st.plotly_chart(fig)
                
                # Detailed positions
                for result in results:
                    with st.expander(f"Positions in {result['sequence_id']}"):
                        if result['positions']:
                            st.write(f"Found at positions: {', '.join(map(str, result['positions']))}")
                        else:
                            st.write("No occurrences found")
    else:
        st.info("Please upload a FASTA file to begin analysis")

def show_sequence_alignment():
    st.header("Sequence Alignment")
    
    if st.session_state.mtdna_objects:
        col1, col2 = st.columns([1, 2])
        
        with col1:
            seq1_id = st.selectbox(
                "Select First Sequence",
                [obj.sequence_id for obj in st.session_state.mtdna_objects],
                key='seq1'
            )
            
            seq2_id = st.selectbox(
                "Select Second Sequence",
                [obj.sequence_id for obj in st.session_state.mtdna_objects],
                key='seq2'
            )
            
            alignment_method = st.radio(
                "Alignment Method",
                ['global', 'local']
            )
        
        if seq1_id != seq2_id:
            seq1 = next(obj for obj in st.session_state.mtdna_objects if obj.sequence_id == seq1_id)
            seq2 = next(obj for obj in st.session_state.mtdna_objects if obj.sequence_id == seq2_id)
            
            aligner = BioPythonAligner()
            
            try:
                # Perform alignment
                alignment_result = aligner.align_sequences(
                    seq1.sequence,
                    seq2.sequence,
                    method=alignment_method
                )
                
                with col2:
                    st.subheader("Alignment Results")
                    
                    # Display alignment statistics
                    st.metric("Alignment Score", f"{alignment_result['score']:.2f}")
                    st.metric("Percent Identity", f"{alignment_result['percent_identity']:.2f}%")
                    st.metric("Number of Gaps", alignment_result['gaps'])
                    
                    # Format and display alignment blocks
                    blocks = aligner.format_alignment(alignment_result)
                    st.subheader("Alignment Visualization")
                    
                    # Create custom CSS for alignment display
                    st.markdown("""
                        <style>
                        .alignment-block {
                            font-family: monospace;
                            white-space: pre;
                            padding: 10px;
                            margin: 10px 0;
                            background-color: #f5f5f5;
                            border-radius: 5px;
                        }
                        .position-header {
                            color: #666;
                            font-size: 0.9em;
                            margin-bottom: 5px;
                        }
                        .sequence {
                            color: #1a1a1a;
                        }
                        .match-line {
                            color: #0066cc;
                        }
                        </style>
                    """, unsafe_allow_html=True)

                    # Display blocks with enhanced formatting
                    window_size = st.slider("Alignment window size", 60, 200, 60)
                    start_position = st.slider("Start position", 0, len(blocks) * 60 - window_size, 0)
                    
                    # Calculate which blocks to show based on window
                    start_block = start_position // 60
                    end_block = (start_position + window_size) // 60 + 1
                    blocks_to_show = blocks[start_block:end_block]

                    for block in blocks_to_show:
                        html = f"""
                        <div class="alignment-block">
                            <div class="position-header">Position {block['position']}</div>
                            <div class="sequence">{block['seq1']}</div>
                            <div class="match-line">{block['match']}</div>
                            <div class="sequence">{block['seq2']}</div>
                        </div>
                        """
                        st.markdown(html, unsafe_allow_html=True)
                    
                    if len(blocks) > 5:
                        st.info("Showing first 5 alignment blocks. Scroll to see more.")
                        
            except Exception as e:
                st.error(f"Error performing alignment: {str(e)}")
        else:
            st.warning("Please select different sequences for alignment")
    else:
        st.info("Please upload a FASTA file to begin analysis")

if __name__ == "__main__":
    main()