# mytheme.jl
# Plots theme
using PlotThemes

sheet_args = Dict{Symbol, Any}([
    :fglegend => plot_color(colorant"#225", 0.1),
    :bglegend => plot_color(:white, 0.9),
    :gridcolor => colorant"#155",
    :minorgridcolor => colorant"#225",
    :framestyle => :grid,
    :minorgrid => true,
    :linewidth => 1.2,
    # :markersize => 6,
    :markerstrokewidth => 0,
    :foreground_color_axis =>:gray30, # tick color
    :tick_direction => :out,
    :foreground_color_text => :gray30, # tick labels
    ]);

mute_palette = [ # bright
    colorant"#c67", # rose
    colorant"#328", # indigo
    colorant"#dc7", # sand
    colorant"#173", # green
    colorant"#8ce", # cyan
    colorant"#825", # wine
    colorant"#4a9", # teal
    colorant"#993", # olive
    colorant"#a49", # purple
    colorant"#ddd", # grey
];

# modified from wong_palette
my_palette = [
    RGB(([230, 159,   0] / 255)...), # orange
    RGB(([ 86, 180, 233] / 255)...), # sky blue
    RGB(([  0, 158, 115] / 255)...), # blueish green
    RGB(([213,  94,   0] / 255)...), # vermillion
    RGB(([  0, 114, 178] / 255)...), # blue
    # RGB(([204, 121, 167] / 255)...), # reddish purple
    ];

wong_palette = [
    RGB(([230, 159,   0] / 255)...), # orange
    RGB(([ 86, 180, 233] / 255)...), # sky blue
    RGB(([  0, 158, 115] / 255)...), # blueish green
    RGB(([240, 228,  66] / 255)...), # yellow
    RGB(([  0, 114, 178] / 255)...), # blue
    RGB(([213,  94,   0] / 255)...), # vermillion
    RGB(([204, 121, 167] / 255)...), # reddish purple
    ];

add_theme(:mytheme, PlotTheme(merge!(
    Dict{Symbol,Any}([
        :palette => my_palette,
        :colorgradient => cgrad(:viridis).colors
        ]),
    sheet_args)
    )
);