import re
import glob
import os

def extract_function_block(text, start_marker):
    """
    Extract the function call block starting at start_marker
    and continuing until the parentheses balance.
    """
    start_index = text.find(start_marker)
    if start_index == -1:
        return None
    
    index = start_index + len(start_marker)
    paren_count = 1  # for the opening parenthesis after the marker
    block = ""
    
    # Loop until we've closed all parentheses
    while index < len(text) and paren_count > 0:
        char = text[index]
        block += char
        if char == '(':
            paren_count += 1
        elif char == ')':
            paren_count -= 1
        index += 1

    return block

def extract_argument_value(block, arg_name):
    """
    Extracts the value of a given argument (arg_name) from a block of text.
    If the value starts with c(, it will read until the matching closing parenthesis.
    Otherwise, it reads until a comma or newline is encountered.
    Additionally, for non-vector values, it trims a trailing ')' if present.
    """
    # Find the start of the argument assignment.
    pattern = re.compile(r'{0}\s*=\s*'.format(re.escape(arg_name)))
    m = pattern.search(block)
    if not m:
        return "NotFound"
    start_index = m.end()
    
    # Skip any whitespace.
    while start_index < len(block) and block[start_index].isspace():
        start_index += 1

    # If the value starts with c( then extract everything until matching ')'.
    if block[start_index:start_index+2] == "c(":
        index = start_index + 2  # start right after "c("
        paren_count = 1
        while index < len(block) and paren_count > 0:
            if block[index] == '(':
                paren_count += 1
            elif block[index] == ')':
                paren_count -= 1
            index += 1
        value = block[start_index:index].strip()
        return value
    else:
        # For simple values, read until a comma or newline is encountered.
        index = start_index
        while index < len(block) and block[index] not in [',', '\n']:
            index += 1
        value = block[start_index:index].strip()
        # If the value ends with a closing parenthesis (e.g. last parameter), remove it.
        if value.endswith(")") and not value.startswith("c("):
            value = value[:-1].strip()
        return value

def main():
    start_marker = "clin_comp = compClinvar_single_algorithm("
    
    # Determine the directory with the R scripts.
    # This script is in Scripts/Survival analysis/ so the R scripts are in:
    rscripts_dir = os.path.join(os.path.dirname(__file__), "..", "single_algorithm")
    
    # Define the output file path.
    # Resources is assumed to be at the same level as Scripts.
    output_filepath = os.path.join(os.path.dirname(__file__), "..", "..", "Resources", "extracted_parameters.txt")
    
    # Write header to the output file.
    header = ("Algorithm\tpdf_level_col_width\tpdf_count_col_width\t"
              "pdf_pval_col_width\tpdf_test_col_width\tpdf_tab_font_size")
    with open(output_filepath, "w", encoding="utf-8") as out_file:
        out_file.write(header + "\n")
    
    # Build a pattern for .R files in the given directory (non-recursive).
    rscript_pattern = os.path.join(rscripts_dir, "*.R")
    for filepath in glob.glob(rscript_pattern):
        # Use utf-8 decoding and replace problematic characters.
        with open(filepath, "r", encoding="utf-8", errors="replace") as file:
            text = file.read()
        
        # Extract the whole function call block.
        block = extract_function_block(text, start_marker)
        if block is None:
            continue

        # Extract desired argument values.
        pdf_level = extract_argument_value(block, "pdf_level_col_width")
        pdf_count = extract_argument_value(block, "pdf_count_col_width")
        pdf_pval  = extract_argument_value(block, "pdf_pval_col_width")
        pdf_test  = extract_argument_value(block, "pdf_test_col_width")
        pdf_font  = extract_argument_value(block, "pdf_tab_font_size")
        
        # Get the "Algorithm" name from the R file name (sans extension).
        algorithm_name = os.path.splitext(os.path.basename(filepath))[0]
        
        # Prepare a tab-separated line.
        line = f"{algorithm_name}\t{pdf_level}\t{pdf_count}\t{pdf_pval}\t{pdf_test}\t{pdf_font}"
        
        with open(output_filepath, "a", encoding="utf-8") as out_file:
            out_file.write(line + "\n")

if __name__ == "__main__":
    main()