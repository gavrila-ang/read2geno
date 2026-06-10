def run_substeps(substep_list, substeps_dict, step_value, check_emoji):
 
    for i, step_name in enumerate(substep_list, start=1):
        
        step_data = substeps_dict[step_name]
        func = step_data['func']
        params = step_data['params']
        desc = step_data['description']

        print(f"--- Starting Step {step_value}.{i}: {step_name} ---")
        print(f"Description: {desc}")
        start = time()

        try:
            func(**params)  # Call the function with unpacked parameters
        except Exception as e:
            print(f"!!! Step {step_value}.{i} ({step_name}) failed with error:\n{e}")
            break
        
        end = time()

        # Calculate time in hours and minutes
        elapsed = end - start
        hours = int(elapsed // 3600)
        minutes = int((elapsed % 3600) // 60)
        print(f"{check_emoji} --- Finished Step {step_value}.{i}: {step_name}. Time taken: {hours}h {minutes}m ---\n")
