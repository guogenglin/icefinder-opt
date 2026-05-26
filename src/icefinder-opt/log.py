# -*- coding: utf-8 -*-
"""
Created on Sun May 17 14:03:07 2026

@author: Genglin Guo
@e-mail: 2019207025@njau.edu.cn
"""

import sys
from inspect import stack
from datetime import datetime


_LOGO = r"""
  ___ ____ _____ _____ ___ _   _ ____  _____ ____  
 |_ _/ ___| ____|  ___|_ _| \ | |  _ \| ____|  _ \ 
  | | |   |  _| | |_   | ||  \| | | | |  _| | |_) |
  | | |___| |___|  _|  | || |\  | |_| | |___|  _ < 
 |___\____|_____|_|   |___|_| \_|____/|_____|_| \_\
"""

LOGO_WIDTH = 50

def generate_gradient_logo(logo_text: str, start_rgb: tuple, end_rgb: tuple) -> str:
    '''
    Generate a text-based logo with a horizontal color gradient from start_rgb to end_rgb.
    '''
    lines = logo_text.strip("\n").split("\n")
    max_len = max(len(line) for line in lines)
    
    colored_logo = []
    
    for y, line in enumerate(lines):
        colored_line = ""
        for x, char in enumerate(line):
            if char == " ":
                colored_line += char
                continue
            
            # Calculate the ratio of the current position to the maximum length for gradient interpolation
            ratio = x / max_len if max_len > 0 else 0
            
            # Interpolate RGB values based on the ratio
            r = int(start_rgb[0] + (end_rgb[0] - start_rgb[0]) * ratio)
            g = int(start_rgb[1] + (end_rgb[1] - start_rgb[1]) * ratio)
            b = int(start_rgb[2] + (end_rgb[2] - start_rgb[2]) * ratio)
            
            # 
            colored_line += f"\033[38;2;{r};{g};{b}m{char}\033[0m"
        colored_logo.append(colored_line)
    
    return "\n".join(colored_logo)

def bold(text: str) -> str:
    '''
    Apply ANSI escape codes to make text bold in the terminal.
    '''
    return f"\033[1m{text}\033[0m"

# Define some colors
CYAN = (0, 255, 255)
MAGENTA = (255, 0, 255)
GOLDEN = (255, 215, 0)

# Generate the text-based logo with gradient
gradient_logo = generate_gradient_logo(_LOGO, start_rgb=CYAN, end_rgb=MAGENTA)

# Define the description message and format it to be centered and bold
description_msg = "MCSMRT: cluster PacBio reads into zOTU/ASV clusters then classify each zOTU/ASV representative."

# Combine the logo and description into a formatted string
formatted_description = f"{gradient_logo}\n\n{bold(description_msg.center(LOGO_WIDTH))}\n"


def log(message: str = '', verbose: bool = True, rjust: int = 20, stack_depth: int = 1):
    '''
    Simple function for logging messages to stderr. Only runs if verbose == True.
    Stack depth can be increased if the parent function name needs to be exposed.
    '''
    if verbose:  # Only build log if verbosity is requested; simple way of controlling log
        sys.stderr.write(f"{datetime.now():%Y-%m-%d %H:%M:%S} {stack()[stack_depth].function:>{rjust}}] {message}\n")


class ICEfinderError(RuntimeError):
    pass