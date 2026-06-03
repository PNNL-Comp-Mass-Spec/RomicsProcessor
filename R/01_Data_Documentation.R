#' ROP Color Palette
#'
#' @description A predefined color palette used throughout RomicsProcessor for consistent visualization
#' of factors and groupings in plots and heatmaps.
#'
#' @format A named character vector of color hex codes
#'
#' @details The ROP_colors palette is automatically loaded when RomicsProcessor is attached and
#' is used as the default color scheme in most plotting functions. Colors are assigned to factor
#' levels in the order they appear in the palette.
#'
#' @examples
#' \dontrun{
#'   # Access the color palette
#'   head(ROP_colors)
#'
#'   # Use in a plot
#'   plot(1:10, col = ROP_colors[1:10], pch = 16, cex = 3)
#' }
#'
"ROP_colors"
