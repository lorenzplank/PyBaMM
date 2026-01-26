import pybamm
import numpy as np
import os
import sys
import importlib.util

def load_parameter_values(path):
    spec = importlib.util.spec_from_file_location("mod", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return pybamm.ParameterValues(mod.get_parameter_values())

def test_summary_variables(parameter_set_name, use_degradation=False):
    print(f"\nTesting parameter set: {parameter_set_name} (Degradation: {use_degradation})")
    try:
        if use_degradation:
            options = {
                "SEI": "ec reaction limited",
                "lithium plating": "reversible",
                "particle mechanics": "swelling and cracking",
                "loss of active material": "stress and reaction-driven",
            }
        else:
            options = {}
            
        model = pybamm.lithium_ion.DFN(options=options)
        
        if parameter_set_name == "Chen2020_mod":
            path = "/Users/lorenzplank/Documents/Projekte/PyBaMM/src/pybamm/input/parameters/lithium_ion/Chen2020_mod.py"
            parameter_values = load_parameter_values(path)
        else:
            parameter_values = pybamm.ParameterValues(parameter_set_name)
        
        # Define a simple experiment
        # Using a shorter experiment to ensure it completes
        experiment = pybamm.Experiment(
            [
                ("Discharge at 1C until 3.5V", "Rest for 5 minutes", "Charge at 1C until 4.1V")
            ]
        )
        
        sim = pybamm.Simulation(model, parameter_values=parameter_values, experiment=experiment)
        sol = sim.solve()
        
        print(f"Simulation successful for {parameter_set_name}")
        print("Summary variables:")
        # Try to access summary variables
        summary = sol.summary_variables
        # summary is a SummaryVariables object
        all_vars = summary.all_variables
        print(f"Found {len(all_vars)} summary variables")
        for key in all_vars[:20]:
            try:
                val = summary[key]
                print(f"  {key}: {val}")
            except Exception as e:
                print(f"  {key}: FAILED ({e})")
            
    except Exception as e:
        print(f"Error for {parameter_set_name}: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    # Test Chen2020 without degradation (likely what works for user)
    test_summary_variables("Chen2020", use_degradation=False)
    # Test Chen2020_mod with degradation (likely what fails for user)
    test_summary_variables("Chen2020_mod", use_degradation=True)
