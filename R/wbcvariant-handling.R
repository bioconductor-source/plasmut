#' Strip a column of a dataframe into a vector
#' @import tidyverse
#' @param data a data frame.
#' @param col.idx the column name or column number
#' @return a vector of the selected column
strip.column <- function(data, col.idx=1){

    x <- data %>% as.list()

    return(x[[col.idx]])
}

#' handle mutations with a non-zero white blood cell mutant read count
#' @import tidyverse
#' @import qpdf
#' @param pdmr vector of distinct mutant reads in plasma for different mutations
#' @param pdr vector of distinct reads in plasma for different mutations
#' @param wdmr vector of distinct mutant reads in wbc for different mutations
#' @param wdr vector of distinct reads in wbc for different mutations
#' @param thresh user specificed threshold for wbcmaf to cluster wbc variant types instead of kmeans clustering
#' @param data data frame with the pdmr, pdr, wdmr, and wdr information as column and rows correspond to different mutations
#' @param cols vector of column names to use (cols[1] refers to pdmr, etc.)
#' @return the vector representing the clusters to which the white blood cell variants belong to
#' @export

stratify_wbc_vars <- function(pdmr, pdr, wdmr, wdr, thresh=NULL, data=NULL, cols=NULL){

    if (!is.null(data)){

        if (!is.null(cols)){

            pdmr <- strip.column(data, cols[1])
            pdr <- strip.column(data, cols[2])
            wdmr <- strip.column(data, cols[3])
            wdr <- strip.column(data, cols[4])

        } else{
            pdmr <- strip.column(data, 1)
            pdr <- strip.column(data, 2)
            wdmr <- strip.column(data, 3)
            wdr <- strip.column(data, 4)
        }
    }

    df.orig <- data.frame(pdmr=pdmr, pdr=pdr, wdmr=wdmr, wdr=wdr) %>% mutate(index=1:n())

    df <- df.orig %>% mutate(wbcmaf=wdmr/wdr, cfdnamaf=pdmr/pdr) %>% mutate(wbcmaf = ifelse(is.nan(wbcmaf), 0, wbcmaf)) %>% filter(wbcmaf > 0)

    if (!is.null(thresh)){

        cluster <- df %>% mutate(mutype = ifelse(wbcmaf >= thresh, 2, 1)) %>% pull(mutype)

    } else{

        k <- df %>% select(wbcmaf) %>% kmeans(., center=2)

        cluster <- k$cluster

    }

    df$wbc_var_type <- cluster

    res <- left_join(df.orig, df %>% select(index, wbc_var_type), by=c("index"="index")) %>% select(-index)

    return(res$wbc_var_type)

}

