def search_and_get_context(file_path, search_string, rows_before=5):
    """
    Search for a string in a file and return the specified number of rows before it.
    
    Args:
        file_path (str): Path to the file to search
        search_string (str): String to search for
        rows_before (int): Number of rows to return before the match
    
    Returns:
        list: List of lines before the match, or None if not found
    """
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        
        # Search for the string
        for i, line in enumerate(lines):
            if search_string in line:
                # Calculate start index for context
                start_idx = max(0, i - rows_before)
                
                # Get the context lines
                context_lines = lines[start_idx:i]
                
                print(f"Found '{search_string}' at line {i + 1}")
                print(f"Showing {len(context_lines)} lines before the match:")
                print("-" * 50)
                
                for j, context_line in enumerate(context_lines):
                    line_num = start_idx + j + 1
                    print(f"{line_num:4d}: {context_line.rstrip()}")
                
                print(f"{i + 1:4d}: {line.rstrip()} <- MATCH")
                print("-" * 50)
                
                return context_lines
        
        print(f"String '{search_string}' not found in file")
        return None
        
    except FileNotFoundError:
        print(f"File '{file_path}' not found")
        return None
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

# Call the function with your specific parameters
if __name__ == "__main__":
    search_and_get_context("/home/htsofie/Desktop/phospho_root/data/blast_dbs/UniprotKB_rat_paper.fasta", "SPYGSRSPFEHSA", rows_before=50)
    