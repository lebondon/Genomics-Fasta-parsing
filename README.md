# Mitochondrial DNA Analysis Tool

A comprehensive web-based application for analyzing mitochondrial DNA sequences, built with Python and Streamlit.

## Overview

This project provides a robust platform for analyzing mitochondrial DNA (mtDNA) sequences from multiple species. It offers various analytical capabilities including:

- FASTA file parsing and visualization
- Sequence analysis (GC content, base composition)
- Motif identification and analysis
- Pairwise sequence alignment
- Interactive visualizations

## Features

- **Data Upload and Parsing**
  - Support for FASTA format files
  - Automatic sequence type detection (DNA/protein)
  - Efficient data organization using Pandas DataFrames

- **Sequence Analysis**
  - Basic statistics (sequence length, GC content)
  - Base composition analysis with visual representations
  - Subsequence extraction
  - Interactive sequence viewer

- **Motif Analysis**
  - Search for specific sequence motifs
  - Visualization of motif distributions
  - Cross-sequence motif comparison

- **Sequence Alignment**
  - Pairwise sequence alignment using BioPython
  - Support for both global and local alignments
  - Visual alignment representation
  - Alignment statistics (score, identity percentage, gaps)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/mtdna-analysis-tool.git
cd mtdna-analysis-tool
```

2. Create and activate a virtual environment (optional but recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install required packages:
```bash
pip install -r requirements.txt
```

## Usage

1. Start the Streamlit application:
```bash
streamlit run app.py
```

2. Open your web browser and navigate to the provided local URL (typically http://localhost:8501)

3. Upload your FASTA file containing mtDNA sequences

4. Navigate through different analysis options using the sidebar:
   - Upload & Overview
   - Sequence Analysis
   - Motif Analysis
   - Sequence Alignment

## Project Structure

```
mtdna-analysis-tool/
├── app.py                 # Main Streamlit application
├── fasta_parser.py        # FASTA file parsing functionality
├── gene_classes.py        # Core classes for sequence analysis
├── requirements.txt       # Project dependencies
└── README.md             # Project documentation
```

## Class Structure

- `Gene`: Base class for genetic sequences
  - `DNASequence`: Class for DNA sequence analysis
  - `AmminoacidsSequence`: Class for protein sequence analysis
  - `BioPythonAligner`: Class for sequence alignment operations

- `FastaParser`: Class for parsing and processing FASTA files

## Development

The project is designed with modularity and extensibility in mind, following object-oriented principles:

- **Inheritance**: Shared behaviors between different sequence types
- **Encapsulation**: Protected sequence data and analysis methods
- **Abstraction**: High-level representation of genomic sequences
- **Polymorphism**: Generic methods for sequence analysis

## Requirements

- Python 3.8+
- Streamlit
- Pandas
- BioPython
- Plotly

## Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Project specification by Andrea Giovanni Nuzzolese
- BioPython library for sequence alignment
- Streamlit for the web interface
