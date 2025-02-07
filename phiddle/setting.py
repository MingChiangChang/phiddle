from dataclasses import dataclass

# Allow fixing of specific setting memeber
# Allow axis sharing (Same xmin xmax for 2 settings)
# Maybe implement with callback
@dataclass
class PlotSettings():
    title: str  = ""
    xlabel: str = "" 
    ylabel: str = ""
    xmin: float = 0.0 
    xmax: float = 1.0 
    ymin: float = 0.0
    ymax: float = 1.0

    @property
    def xlims(self):
        return (self.xmin, self.xmax)

    @property
    def ylims(self):
        return (self.ymin, self.ymax)


@dataclass
class HeatmapPlotSettings(PlotSettings): 
    vmin: float = 0.0
    vmin: float = 1.0
    cmap: str = "" # FIXME: Find a defual cmap


@dataclass
class ScatterPlotSettings(PlotSettings): 
    markersize: float = 1.0
    cmap: str = ""

    
def set_plot(ax, settings: PlotSettings, callback=None):
    ax.set_title(settings.title)
    ax.set_xlabel(settings.xlabel)
    ax.set_ylabel(settings.ylabel)
    ax.set_xlim(settings.xlims)
    ax.set_ylim(settings.ylims)

    if callback is not None:
        callback()
