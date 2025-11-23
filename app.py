from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
import subprocess
import os
import json
from pathlib import Path

app = Flask(__name__)
CORS(app)  # Enable CORS for local React app

# Configuration
MATLAB_PATH = r"C:\Program Files\MATLAB\R2025a\bin\matlab.exe"  # Update this to your MATLAB path
CIRCUIT_ANALYSIS_SCRIPT = "Circuit_Analysis.m"
NETLIST_FILE = "circuit_netlist.txt"
RESULTS_FILE = "Results.txt"


@app.route('/api/health', methods=['GET'])
def health_check():
    """Check if backend is running"""
    return jsonify({"status": "online", "message": "Backend is running"})


@app.route('/api/simulate', methods=['POST'])
def simulate_circuit():
    """Receive netlist from frontend and run MATLAB simulation"""
    try:
        data = request.json
        netlist = data.get('netlist', [])
        tf = data.get('tf', 1.0)  # Final time for transient analysis

        if not netlist:
            return jsonify({"error": "No netlist provided"}), 400

        # Write netlist to file
        with open(NETLIST_FILE, 'w') as f:
            for line in netlist:
                f.write(' '.join(str(x) for x in line) + '\n')

        # Create a MATLAB script that runs the simulation
        matlab_command = f"""
        cd '{os.getcwd()}';
        tf = {tf};
        Circuit_Analysis;
        exit;
        """

        # Save MATLAB command to temporary file
        with open('run_simulation.m', 'w') as f:
            f.write(matlab_command)

        # Run MATLAB in batch mode
        print(f"Running MATLAB simulation with netlist: {NETLIST_FILE}")

        # Method 1: Using matlab command line (Windows)
        process = subprocess.Popen(
            [MATLAB_PATH, '-batch', f"run('{os.path.abspath('run_simulation.m')}')",'-nosplash'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        stdout, stderr = process.communicate(timeout=180)

        # Check if simulation completed
        if not os.path.exists(RESULTS_FILE):
            return jsonify({
                "error": "Simulation failed - no results generated",
                "matlab_output": stdout,
                "matlab_error": stderr
            }), 500

        # Read results
        with open(RESULTS_FILE, 'r') as f:
            results = f.read()

        # Check for generated figures
        figures = []
        # --- CHANGE: Look for .png files instead of .fig ---
        figure_files = ['Voltages_graph.png', 'Currents_graph.png', 'Cap_curr.png', 'Res_C.png']
        for fig_file in figure_files:
            if os.path.exists(fig_file):
                figures.append(fig_file)

        # Clean up temporary file
        if os.path.exists('run_simulation.m'):
            os.remove('run_simulation.m')

        return jsonify({
            "success": True,
            "results": results,
            "figures": figures,
            "matlab_output": stdout
        })

    except subprocess.TimeoutExpired:
        return jsonify({"error": "MATLAB simulation timed out"}), 500
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route('/api/simulate-python', methods=['POST'])
def simulate_circuit_python():
    """Alternative: Run MATLAB using Python's MATLAB Engine (if installed)"""
    try:
        import matlab.engine

        data = request.json
        netlist = data.get('netlist', [])

        if not netlist:
            return jsonify({"error": "No netlist provided"}), 400

        # Write netlist to file
        with open(NETLIST_FILE, 'w') as f:
            for line in netlist:
                f.write(' '.join(str(x) for x in line) + '\n')

        # Start MATLAB engine
        print("Starting MATLAB engine...")
        eng = matlab.engine.start_matlab()

        # Run the circuit analysis
        print("Running Circuit_Analysis.m...")
        eng.eval(f"Circuit_Analysis", nargout=0)

        # Stop MATLAB engine
        eng.quit()

        # Read results
        if not os.path.exists(RESULTS_FILE):
            return jsonify({"error": "Simulation failed - no results generated"}), 500

        with open(RESULTS_FILE, 'r') as f:
            results = f.read()

        # Check for generated figures
        figures = []
        # --- CHANGE: Look for .png files instead of .fig ---
        figure_files = ['Voltages_graph.png', 'Currents_graph.png', 'Cap_curr.png', 'Res_C.png']
        for fig_file in figure_files:
            if os.path.exists(fig_file):
                figures.append(fig_file)

        return jsonify({
            "success": True,
            "results": results,
            "figures": figures
        })

    except ImportError:
        return jsonify({
            "error": "MATLAB Engine for Python not installed. Use /api/simulate endpoint instead."
        }), 500
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route('/api/results', methods=['GET'])
def get_results():
    """Get the latest simulation results"""
    try:
        if not os.path.exists(RESULTS_FILE):
            return jsonify({"error": "No results available"}), 404

        with open(RESULTS_FILE, 'r') as f:
            results = f.read()

        return jsonify({"results": results})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route('/api/download/<filename>', methods=['GET'])
def download_file(filename):
    """Download generated files (results, figures)"""
    try:
        if os.path.exists(filename):
            return send_file(filename, as_attachment=True)
        else:
            return jsonify({"error": "File not found"}), 404
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route('/api/config', methods=['GET'])
def get_config():
    """Get current configuration"""
    return jsonify({
        "matlab_path": MATLAB_PATH,
        "matlab_installed": os.path.exists(MATLAB_PATH),
        "circuit_analysis_exists": os.path.exists(CIRCUIT_ANALYSIS_SCRIPT)
    })


@app.route('/api/config', methods=['POST'])
def update_config():
    """Update MATLAB path"""
    global MATLAB_PATH
    data = request.json
    new_path = data.get('matlab_path')

    if new_path and os.path.exists(new_path):
        MATLAB_PATH = new_path
        return jsonify({"success": True, "matlab_path": MATLAB_PATH})
    else:
        return jsonify({"error": "Invalid MATLAB path"}), 400

# --- NEW FUNCTION ---
@app.route('/api/figures/<filename>')
def get_figure(filename):
    """Serve a generated figure file"""
    try:
        # Basic security: only allow .png files from the current directory
        if filename.endswith('.png') and os.path.exists(filename):
            # send_file without 'as_attachment=True' displays the image
            return send_file(filename, mimetype='image/png')
        else:
            return jsonify({"error": "File not found or not allowed"}), 404
    except Exception as e:
        return jsonify({"error": str(e)}), 500
# --------------------

if __name__ == '__main__':
    print("=" * 60)
    print("GSpice Circuit Simulator - Backend Server")
    print("=" * 60)
    print(f"MATLAB Path: {MATLAB_PATH}")
    print(f"MATLAB Installed: {os.path.exists(MATLAB_PATH)}")
    print(f"Circuit Analysis Script: {os.path.exists(CIRCUIT_ANALYSIS_SCRIPT)}")
    print("=" * 60)
    print("\nServer starting on http://localhost:5000")
    print("Press Ctrl+C to stop the server\n")

    app.run(debug=True, port=5000, host='0.0.0.0')