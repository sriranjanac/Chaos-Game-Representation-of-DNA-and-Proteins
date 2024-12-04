# Import necessary libraries
from flask import Flask, render_template, request, redirect, url_for  
import os 
import matplotlib.pyplot as plt  
from werkzeug.utils import secure_filename 
from Bio import SeqIO  
import numpy as np  

# Initialize the Flask app and configure the upload folder
app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = os.path.join('static', 'uploads')  # Path to save uploaded files
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)  # Create upload directory if it doesn't exist

# Processing the fasta files
def process_fasta(file_path):
    """A function to read and process the uploaded FASTA file."""
    with open(file_path, 'r') as file:
        content = file.read()
    return content

# CGR Function for DNA sequence
def generate_dna_cgr(sequence, output_path):
    # Define coordinates for each DNA base
    base_mapping = {
        'A': (0, 0),
        'C': (1, 0),
        'G': (1, 1),
        'T': (0, 1)
    }

    x, y = 0.5, 0.5  # Start from the center of the plot
    x_coords, y_coords = [x], [y]  # Lists to store coordinates

    # Iterate over each base in the DNA sequence and update coordinates
    for base in sequence:
        if base in base_mapping:
            x = (x + base_mapping[base][0]) / 2
            y = (y + base_mapping[base][1]) / 2
            x_coords.append(x)
            y_coords.append(y)

    # Plot the CGR for DNA
    plt.figure(figsize=(5, 5))
    plt.scatter(x_coords, y_coords, s=0.1, color='black')  # Plot each point with a small size

    # Draw boundaries around the plot
    line_thickness = 0.6
    plt.plot([0, 1], [0, 0], color='black', linewidth=line_thickness)  # Bottom boundary
    plt.plot([0, 1], [1, 1], color='black', linewidth=line_thickness)  # Top boundary
    plt.plot([0, 0], [0, 1], color='black', linewidth=line_thickness)  # Left boundary
    plt.plot([1, 1], [0, 1], color='black', linewidth=line_thickness)  # Right boundary

    # Add labels for each corner
    plt.annotate('A', (0, 0), textcoords="offset points", xytext=(-20, -20), fontsize=12, color='red', fontweight="bold")
    plt.annotate('C', (1, 0), textcoords="offset points", xytext=(5, -20), fontsize=12, color='red', fontweight="bold")
    plt.annotate('G', (1, 1), textcoords="offset points", xytext=(0, 15), fontsize=12, color='red', fontweight="bold")
    plt.annotate('T', (0, 1), textcoords="offset points", xytext=(-15, 15), fontsize=12, color='red', fontweight="bold")

    plt.title("CGR of DNA", fontweight="bold", fontsize=12)
    plt.axis('off')
    plt.savefig(output_path)  # Save the plot to the specified output path
    plt.close()  # Close the plot to free up memory

# CGR Function for Protein sequence
def generate_protein_cgr(sequence, output_path):
    # Generate coordinates for a dodecagon for 12 groups of amino acids
    angles = np.linspace(0, 2 * np.pi, 13)[:-1]  # 12 points equally spaced around a circle
    angles = angles[::-1]  # Reverse order for clockwise layout
    radius = 0.5
    dodecagon_coords = [(radius * np.cos(angle) + 0.5, radius * np.sin(angle) + 0.5) for angle in angles]

    # Map each group of amino acids to a vertex of the dodecagon
    group_mapping = {
        'ILVM': dodecagon_coords[8], 'RK': dodecagon_coords[9], 'DE': dodecagon_coords[10],
        'N': dodecagon_coords[11], 'Q': dodecagon_coords[0], 'H': dodecagon_coords[1],
        'ST': dodecagon_coords[2], 'P': dodecagon_coords[3], 'AG': dodecagon_coords[4],
        'C': dodecagon_coords[5], 'FY': dodecagon_coords[6], 'W': dodecagon_coords[7]
    }

    x, y = 0.5, 0.5  # Start from the center
    x_coords, y_coords = [x], [y]  # Lists to store coordinates

    # Map each amino acid in the sequence to its corresponding group and update coordinates
    for amino_acid in sequence:
        for group, coord in group_mapping.items():
            if amino_acid in group:
                x = (x + coord[0]) / 2
                y = (y + coord[1]) / 2
                x_coords.append(x)
                y_coords.append(y)
                break

    # Plot the CGR for proteins
    plt.figure(figsize=(6, 6))
    plt.scatter(x_coords, y_coords, s=0.1, color='black')  # Plot each point with a small size

    # Connect the vertices of the dodecagon with lines
    for i, (coord, next_coord) in enumerate(zip(dodecagon_coords, dodecagon_coords[1:] + dodecagon_coords[:1])):
        plt.plot([coord[0], next_coord[0]], [coord[1], next_coord[1]], 'k-', lw=0.6)

    # Add labels for each group at the vertices
    group_labels = ["Q", "H", "ST", "P", "AG", "C", "FY", "W", "ILVM", "RK", "DE", "N"]
    label_offsets = [(10, -10), (0, -20), (-10, -20), (0, -20), (-15, -20), (-15, -10), 
                     (-10, 10), (-5, 8), (-5, 5), (10, -5), (10, -10), (10, -10)]

    for i, (coord, label, offset) in enumerate(zip(dodecagon_coords, group_labels, label_offsets)):
        plt.annotate(label, coord, textcoords="offset points", xytext=offset, fontsize=12, color='red', fontweight="bold")

    plt.title("CGR of Proteins", fontweight="bold", fontsize=12)
    plt.axis('off')
    plt.savefig(output_path)  # Save the plot to the specified output path
    plt.close()  # Close the plot to free up memory

# Home route to display the upload form
@app.route('/')
def home():
    return render_template('index.html')  # Render the home page template

# Route to process the uploaded file and display the result
@app.route('/result', methods=['POST'])
def result():
    if 'fasta_file' not in request.files:
        return redirect(request.url)  # Redirect if no file is uploaded

    fasta_file = request.files['fasta_file']
    if fasta_file.filename == '':
        return redirect(request.url)  # Redirect if the filename is empty

    if fasta_file:
        filename = secure_filename(fasta_file.filename)  # Securely handle filename
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        fasta_file.save(filepath)  # Save the uploaded file

        # Parse and read the sequence from the FASTA file
        record = next(SeqIO.parse(filepath, "fasta"))
        sequence = str(record.seq).upper()

        # Determine sequence type and generate CGR
        sequence_type = request.form['sequence_type']
        output_image = os.path.join(app.config['UPLOAD_FOLDER'], 'cgr.png')

        if sequence_type == 'DNA':
            generate_dna_cgr(sequence, output_image)  # Generate CGR for DNA
        elif sequence_type == 'Protein':
            generate_protein_cgr(sequence, output_image)  # Generate CGR for Protein

    return render_template('result.html', filename=filename, sequence_type=sequence_type, cgr_image='uploads/cgr.png')  # Render the result page

# Run the app if the script is executed directly
if __name__ == '__main__':
    app.run(debug=True)
