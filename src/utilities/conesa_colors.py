from matplotlib.colors import LinearSegmentedColormap, to_hex

CONESA_COLORS = {"green": "#15918A", "orange": "#F58A53", "yellow": "#FDC659", "blue": "#74CDF0",
                 "purple": "#9F7BB8", "pink": "#FDA3D1", "magenta": "#EE446F"}
CONESA_PALETTES = {name: [CONESA_COLORS[c] for c in cols] for name, cols in {
    "main": ["green", "orange", "yellow"],
    "nature": ["green", "yellow"],
    "sunshine": ["yellow", "orange"],
    "hot": ["orange", "magenta"],
    "warm": ["blue", "magenta"],
    "cold": ["blue", "purple"],
    "complete": ["green", "orange", "yellow", "blue", "pink", "purple", "magenta"]}.items()}


def conesa_cmap(palette="main", reverse=False):
    '''Continuous colormap interpolating the colours of a Conesa palette, like conesa_pal in R'''
    colors = CONESA_PALETTES[palette][::-1] if reverse else CONESA_PALETTES[palette]
    return LinearSegmentedColormap.from_list(f"conesa_{palette}", colors)


def color_conesa(n, palette="complete", reverse=False):
    '''n colours of a Conesa palette, like colorConesa in R: the first n of the palette when it has
    enough, otherwise n evenly spaced colours of its interpolation (both ends included)'''
    colors = CONESA_PALETTES[palette][::-1] if reverse else CONESA_PALETTES[palette]
    if n <= len(colors):
        return colors[:n]
    cmap = conesa_cmap(palette, reverse)
    return [to_hex(cmap(i / (n - 1))) for i in range(n)]