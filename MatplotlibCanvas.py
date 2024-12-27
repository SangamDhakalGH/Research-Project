from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class MatplotlibCanvas(FigureCanvas):
    """ This small class is used to set up a Matplotlib figure and axes. """
    
    def __init__(self, parent=None, width = 10, height = 10, dpi = 100):
        self.fig  = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MatplotlibCanvas, self).__init__(self.fig)
        
    def refresh(self) -> None:
        """ Update the graph to include any new data that was applied to it. """
        
        self.draw()
    
    def wipe(self) -> None:
        """ Clear all plots and labels from the graph. """
        
        self.axes.cla()
        
    def clear_lines(self) -> None:
        """ Remove the plots from the graph, keeping the labels and zoom. """
        
        for line in self.fig.gca().lines:
            line.remove()
        self.refresh()

    def draw_plot(self, x=[], y=[], color='tab:blue'):
        """
        Plot data to the graph.

        # Inputs
        - x (any): The x-axis data.
        - y (any): The y-axis data.
        - color (str, optional): Color of the plot. Defaults to 'tab:blue'.
        
        # Returns
        - line (matplotlib.lines.Line2D): A list containing the lines of the plot.
        """
        
        if list(x) == []: line = self.axes.plot(y, color=color,zorder=2)
        else:       line = self.axes.plot(x, y, color=color,zorder=2)
        self.refresh()
        
        return line
        
    def draw_scatter(self, x, y, color='tab:blue'):
        """
        Produce a scatter plot.
        
        # Inputs
        - x (any): The x-axis data.
        - y (any): The y-axis data.
        - color (str, optional): Color of the plot. Defaults to 'tab:blue'.
            
        # Returns
        - line (matplotlib.lines.Line2D): A list containing the lines of the plot.
        """
        
        line = self.axes.scatter(x, y, color=color)
        self.refresh()
        
        return line
    
    def draw_baseline(self, x, y):
        """
        Draws a horizontal line for the baseline of the data.
        """
        self.axes.plot(x,([y]*len(x)), color='grey', linestyle=':', alpha=0.5, zorder=2.5)
        self.refresh()