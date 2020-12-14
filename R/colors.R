#### Just a bunch of colors ####
many <- c("#8E13A2", "#198749", "#2d62a3", "#ffbb00", "#ED49B1", "#Ff7708", "#573794", "#15889C", "#Ae0031", "#9e9e9e", "#795548", "#CD6155", "#F7DC6F", "#7DCEA0", "#85C1E9", "#EB984E", "#ead9d5", "#03cd4A", "#CDDC39", "#e0594b", "#c76cde", "#24B177", "#8D6E63", "#486ff7", "#6300b5", "#88e200", "#012824", "#0d3290", "#a347fb", "#54fc7a", "#eb1388", "#b0978d", "#fe52cf", "#83f1f6", "#f1f847", "#2b1dfc", "#6c6f15", "#6ca05c", "#7788cd", "#f502f3", "#0dc290", "#fa0e03", "#3caa0a", "#befc8d", "#08f8eb", "#b1cd3f", "#d6a5fa", "#ce606c", "#ab1eba", "#6ecc9f", "#054ddc", "#486ff7", "#854f49", "#f22B21", "#3a0e43", "#225805", "#37d160", "#e4b974", "#a8bade", "#47edd1", "#f47a92", "#c76cde", "#9106eb", "#81aa20", "#d7fdfd", "#5deb2e", "#f82745", "#6435e0", "#027ffe", "#8e3101", "#16f648", "#1c15bc", "#8be46e", "#8d6fa0", "#e68fc6", "#058ca9", "#9e018a", "#bdfd0b", "#b22760", "#2bf49f", "#cb9348", "#9d8303", "#c251a1", "#46adaf", "#a3e3af", "#22bb34", "#6ea3fa", "#260374", "#1c3854", "#405d37", "#c21df3", "#fcea92", "#537f88", "#fd4c18", "#f2d71e", "#fd4c7a")
gyrdpu <- grDevices::colorRampPalette(c('#e5e7e9', RColorBrewer::brewer.pal(n=9, name="RdPu")))(100)
gyrdpu_flat <- grDevices::colorRampPalette(c('#e5e7e9', RColorBrewer::brewer.pal(n=9, name="RdPu")), bias = 0.5)(100)
inferno <- grDevices::colorRampPalette(viridis::viridis(n=10, option = 'inferno'))(100)
inferno_flat <- grDevices::colorRampPalette(viridis::viridis(n=10, option = 'inferno'), bias = 0.5)(100)
blues <- colorRampPalette(c('white', brewer.pal(n=9, name='Blues')))(100)
blues_flat <- colorRampPalette(c('white', brewer.pal(n=9, name='Blues')), bias=0.8)(100)
ylorrd <- colorRampPalette(brewer.pal(n=9, name='YlOrRd'))(100)
ylorrd_flat <- colorRampPalette(brewer.pal(n=9, name='YlOrRd'), bias=0.8)(100)

#### Themes for plots ####
feature_color_scale <- scale_color_gradientn(colours=gyrdpu)
feature_fill_scale <- scale_fill_gradientn(colours=gyrdpu)
feature_theme <- theme(legend.position = 'none')
dr_theme <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank()
)
dr_guides <- guides(
    color = guide_legend(override.aes = list(size = 5, alpha=1)),
    fill = guide_legend(override.aes = list(size = 5, alpha=1)),
    alpha = guide_legend(override.aes = list(size = 5))
)

#### Structure colors ####
struct_names_custom2 <- c(
    'olfactory bulb',
    'pallium',
    'subpallium',
    'preoptic telencephalon',
    'telencephalo-hypothalamic transition area',
    'hypothalamus',
    'diencephalon',
    'midbrain',
    'prepontine hindbrain',
    'pontine hindbrain',
    'promedullary hindbrain',
    'medullary hindbrain',
    'spinal cord',
    'telencephalic vesicle (alar plate)',
    'telencephalic vesicle (roof plate)'
)

struct_colors_custom2 <- c(
    '#f48fb1',
    '#ad1457',
    '#7b1fa2',
    '#5e35b1',
    '#424242',
    '#ba68c8',
    '#9fa8da',
    '#0097a7',
    '#1b5e20',
    '#aed581',
    '#558b2f',
    '#c0ca33',
    '#fbc02d',
    '#757575',
    '#bdbdbd'
)

struct_symbol_custom2 <- c(
    'OB',
    'Pall',
    'SPall',
    'POTel',
    'TelH',
    'Hyp',
    'D',
    'M',
    'PPH',
    'PH',
    'PMH',
    'MH',
    'SpC',
    'TelA',
    'TelR'
)

names(struct_colors_custom2) <- struct_names_custom2

struct_custom2 <- tibble(
    annot = struct_names_custom2,
    colors = struct_colors_custom2,
    symbol = struct_symbol_custom2
)

#### Brainspan colors ####
bs_names <- c('NCx', 'HIP', 'AMY', 'GE', 'STR', 'DTH', 'URL', 'CB')
bs_colors <- c('#ad1457', '#f48fb1', '#884ea0', '#7b1fa2', '#5e35b1',
               '#9fa8da', '#1b5e20', '#aed581')
names(bs_colors) <- bs_names


bs_age_colors <- c(
    '#BDBDBD', '#636363',
    '#E999DA', '#883D79',
    '#E9A0A6', '#E55868', '#8E3637',
    '#EFB935', '#916C22', '#AED15A', '#5D7B2E',
    '#7372D7', '#3A397D')
names(bs_age_colors) <- c(
    8,9,
    12,13,
    16,17,19,
    21,24,25,26,
    35,37
)

#### Mousebrain colors ####
mb_names <- c('Forebrain', 'Midbrain', 'Hindbrain', 'Head')
mb_colors <- c('#ab1673', '#007A92', '#239b56', '#fbc02d')
names(mb_colors) <- mb_names






