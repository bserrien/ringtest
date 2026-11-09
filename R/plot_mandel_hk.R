#' plot_mandel_hk
#'
#' @description
#' Plot Mandel's h and k statistics (ISO 8.3.2)
#'
#' @param mandel outcome of function `mandel_hk`
#'
#' @import ggplot2
#' @importFrom ggpubr ggarrange
#'
#' @export
plot_mandel_hk <- function(mandel) {

  mh <- mandel |>
    ggplot(aes(lab, h, fill = factor(h_sign))) +
    facet_wrap(~ level, nrow = 1) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    geom_text(data = subset(mandel, h_sign == 1), aes(label = lab), size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 3)) +
    labs(y = "Mandel h (outlying mean)", x = "")

  mk <- mandel |>
    ggplot(aes(lab, k, fill = factor(k_sign))) +
    facet_wrap(~ level, nrow = 1) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    geom_text(data = subset(mandel, k_sign == 1), aes(label = lab), size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 3)) +
    labs(y = "Mandel k (outlying variance)", x = "lab")

  ggpubr::ggarrange(mh, mk, nrow = 2)
}
