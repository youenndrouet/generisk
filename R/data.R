#' Eriscam french dataset of 236 MLH1 Lynch Syndrome families
#'
#' Data on 236 Lynch Syndrome (LS) french families.
#' First analysis of the full eriscam dataset comprising 537 LS families with the GRL method is published
#' here (Bonadona et al. JAMA 2011) <https://jamanetwork.com/journals/jama/fullarticle/900645>
#'
#' @format ## `eriscam_mlh1`
#' A data frame with 4,703 rows and 22 columns:
#' \describe{
#'   \item{CENTER_NAME}{Center identification number (de-identified)}
#'   \item{FAMILY_ID}{Family identification number (de-identified)}
#'   \item{PERSON_ID}{Individual identification number (unique within a given family)}
#'   \item{FATHER_ID, MOTHER_ID}{Father's and Mother's identification numbers, unique within a given family. Founders have value (0,0)}
#'   \item{SEX}{1: Male; 2: Female}
#'   \item{PROBAND_FLAG}{1: proband (first person being tested), also called the "Index case"}
#'   \item{MLH1_STATUS}{Genotype of the MLH1 gene; 1: deleterious mutation found; 0: no mutation found ; 4 : unknown mutation status }
#'   \item{MSH2_STATUS}{Genotype of the MSH2 gene; 1: deleterious mutation found; 0: no mutation found ; 4 : unknown mutation status }
#'   \item{MSH6_STATUS}{Genotype of the MSH6 gene; 1: deleterious mutation found; 0: no mutation found ; 4 : unknown mutation status }
#'   \item{AGE_AT_LAST_NEWS}{Age at last news (years)}
#'   \item{COLORECTUM}{Age at colorectal cancer (years)}
#'   \item{ENDOMETRIUM}{Age at endometrial cancer (years)}
#'   \item{STOMACH}{Age at stomach cancer (years)}
#'   \item{SMALL_BOWEL}{Age at small bowel cancer (years)}
#'   \item{OVARY}{Age at ovary cancer (years)}
#'   \item{UROTHELIUM}{Age at urothelial cancer (years)}
#'   \item{BILIARY_TRACT}{Age at biliary tract cancer (years)}
#'   \item{FIRST_COLONOSCOPY}{Age at first colonoscopy (years)}
#'   \item{HYSTERECTOMY}{Age at hysterectomy (years)}
#'   \item{OOPHORECTOMY}{Age at oophorectomy (years)}
#'   \item{TOTAL_COLECTOMY}{Age at total colectomy (years)}
#' }
#' @source <https://jamanetwork.com/journals/jama/fullarticle/900645>
"eriscam_mlh1"
