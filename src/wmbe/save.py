"""
---------------------------------------------------------------------

Tools to export plots to PDF etc format files.

---------------------------------------------------------------------

Requires Python packages/modules:
  -  :mod:`os`

---------------------------------------------------------------------

"""

import os 

def create_plots_dir(results_dir_name='Results1d'):
    """
    Create an output directory for plots.
    
    Throws an exception if the directory cannot be created.
    Returns quietly if the directory already exists.
    
    Attributes:
        results_dir_name (str) : name of directory
    
    Returns:
        str: path to directory (see :mod:`os.path`)

    """
    results_dir = os.path.join('.',results_dir_name)
    try:
        if not os.path.exists(results_dir):
            os.mkdir(results_dir)
        else:
            return results_dir
    except OSError:  
        print('Cannot create results directory')
        raise
    except:  
        raise
    return results_dir


def export_plots(fdict, results_dir, file_type='pdf'):
    """
    Export plots to PDFs or other format files
    
    Attributes:
        fdict (dict) : figures dictionary
        results_dir (str) : name of output directory
        file_type (str) : file format

    """
    for fdi in list(fdict.items()):
        export_plot(*fdi, results_dir, file_type=file_type)


def export_plot(fig_name,fig,results_dir, file_type='pdf'):
    """
    Export plot to PDF or other format file
    
    Attributes:
        fig_name (str) : name to be used for file (extension auto-appended)
        fig (obj) : figure object
        results_dir (str) : name of output directory
        file_type (:mod:`MatPlotLib/Pyplot <matplotlib.pyplot>`) : file format

    """
    fig_name += '.'+file_type.lower()
    try:
        fig.savefig(os.path.join(results_dir,fig_name),
                   bbox_inches = 'tight', pad_inches = 0)
        print('Exported "'+fig_name+'"')
    except OSError:  
        print('Failed to export figure "'+fig_name+'"')
        raise
    except:
        raise