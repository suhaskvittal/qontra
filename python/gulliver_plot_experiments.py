"""
    author: Suhas Vittal
    date:   23 August 2022
"""

from svittal_plotlib import *
from gulliver_data_io import *

DPI = 300

# Setup light theme for figures.
light_theme()

def plot_mwpm_timing_analysis_file(code_dist=3):
    file_path = '../data/mwpm_timing_d=%d.txt' % code_dist
    output_path = '../data/figures/mwpm_timing_d=%d.png' % code_dist

    df = read_timing_analysis_file(file_path, max_time=3000)
    fig = mkhist(df, x_name='Time Taken', y_name='count',\
            category_key='Hamming Weight', palette='dark',\
            binwidth=50, multiple='dodge')
    add_title(fig,\
            'Execution Time Distribution for '
            'Distance %d Surface Code, MWPM' % code_dist)
    plt.savefig(output_path, dpi=DPI) 
    plt.show()

def plot_gulliver_timing_analysis_file(code_dist=3):
    file_path = '../data/bfu_timing_d=%d.txt' % code_dist
    output_path = '../data/figures/bfu_timing_d=%d.png' % code_dist

    df = read_timing_analysis_file(file_path, max_time=3000)
    fig = mkhist(df, x_name='Time Taken', y_name='count',\
            category_key='Hamming Weight', palette='dark',\
            binwidth=10, multiple='dodge')
    add_title(fig,\
            'Execution Time Distribution for '
            'Distance %d Surface Code, Gulliver' % code_dist)
    plt.savefig(output_path, dpi=DPI) 
    plt.show()

if __name__ == '__main__':
    for d in [3, 5, 7, 9]:
        plot_mwpm_timing_analysis_file(code_dist=d)
        plot_gulliver_timing_analysis_file(code_dist=d)

