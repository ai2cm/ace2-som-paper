import matplotlib as mpl
import matplotlib.pyplot as plt


FONTSIZE = 8


def configure_style(fontsize=FONTSIZE, legend_fontsize=FONTSIZE):
    mpl.rcParams["font.size"] = fontsize
    mpl.rcParams["axes.titlesize"] = fontsize
    mpl.rcParams["axes.labelsize"] = fontsize
    mpl.rcParams["xtick.labelsize"] = fontsize
    mpl.rcParams["ytick.labelsize"] = fontsize
    mpl.rcParams["legend.fontsize"] = legend_fontsize
    mpl.rcParams["figure.titlesize"] = fontsize
    mpl.rcParams["savefig.transparent"] = True


def width_constrained_plot(
    nrows,
    ncols,
    width,
    aspect,
    left_pad=0.5,
    right_pad=0.5,
    bottom_pad=0.5,
    top_pad=0.5,
    horizontal_pads=0.2,
    vertical_pads=0.2,
    cbar_pad=0.2,
    cbar_thickness=0.125,
    axes_kwargs={},
):
    fig = plt.figure()

    if isinstance(horizontal_pads, float):
        horizontal_pads = [horizontal_pads for _ in range(ncols - 1)]
    if isinstance(vertical_pads, float):
        vertical_pads = [vertical_pads for _ in range(nrows - 1)]

    panel_width = (width - (left_pad + right_pad + sum(horizontal_pads))) / ncols
    panel_height = aspect * panel_width

    height = (
        nrows * panel_height
        + bottom_pad
        + top_pad
        + sum(vertical_pads)
        + cbar_pad
        + cbar_thickness
    )

    relative_panel_width = panel_width / width
    relative_panel_height = panel_height / height

    axes = []
    for row in range(nrows):
        bottom = (
            row * relative_panel_height
            + (sum(vertical_pads[:row]) + bottom_pad + cbar_pad + cbar_thickness)
            / height
        )
        for col in range(ncols):
            left = (
                col * relative_panel_width
                + (sum(horizontal_pads[:col]) + left_pad) / width
            )
            position = (left, bottom, relative_panel_width, relative_panel_height)
            axes.append(fig.add_axes(position, **axes_kwargs))

    caxes = []
    for col in range(ncols):
        bottom = bottom_pad / height
        left = (
            col * relative_panel_width + (sum(horizontal_pads[:col]) + left_pad) / width
        )
        position = (left, bottom, relative_panel_width, cbar_thickness / height)
        caxes.append(fig.add_axes(position))

    fig.set_size_inches(width, height)
    return fig, axes, caxes
