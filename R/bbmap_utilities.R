#' BBMap Install
#'
#' @param url
#' @param destdir
#'
#' @return
#' @export
#'
#' @examples
bbmap_install <- function(url, destdir = "bin"){

  if(missing(url)){
    url <- ("https://downloads.sourceforge.net/project/bbmap/BBMap_38.57.tar.gz")
  }

  if(!dir.exists(destdir)){dir.create(destdir)} # Create first directory
  if(dir.exists(paste0(destdir,"/bbmap"))){unlink(paste0(destdir,"/bbmap"),recursive=TRUE)} # Remove old version

  destfile <- file.path(destdir, basename(url))
  if(exists(destfile)){file.remove(destfile)} # Remove old zip file

  utils::download.file(url, destfile = destfile)
  utils::untar(destfile,exdir = destdir)  ## check contents
  file.remove(destfile)
}


# Trim primers ------------------------------------------------------------


#' Trim primers using BBDuk
#'
#' @param install
#' @param fwd
#' @param rev
#' @param primers
#' @param restrictleft
#' @param outpath
#' @param ktrim
#' @param ordered
#' @param kmer
#' @param mink
#' @param tpe
#' @param hdist
#' @param copyundefined
#' @param overwrite
#' @param quality
#' @param maxlength
#'
#' @return
#' @export
#'
#' @examples
bbtools_trim <- function(install=NULL, fwd, rev=NULL, primers,
                         restrictleft=NULL,outpath="bbduk",ktrim="l", ordered=TRUE, kmer=NULL,mink=FALSE,tpe=TRUE, hdist=0,copyundefined=TRUE,
                         overwrite=TRUE,quality=TRUE,maxlength=NULL){
  nsamples <- length(fwd)

  bbduk <- function(install=NULL, fwd, rev=NULL, primers,
                    restrictleft=NULL,outpath="bbduk",ktrim="l", ordered=TRUE, kmer=NULL,mink=FALSE,tpe=TRUE, hdist=0,copyundefined=TRUE,
                    overwrite=TRUE,quality=TRUE,maxlength=NULL){


    install= paste0(install,"/current jgi.BBDuk")

    in1=paste0( "in=",fwd)
    if(!is.null(rev)){in2=paste0( "in2=",rev)}  else (in2="")

    if(!is.null(primers)){literal=paste0("literal=",paste0(primers,collapse=",")) } else (stop("Primer sequences are required for trimming"))

    out = paste0("out=",dirname(fwd),"/",outpath,"/", basename(fwd) ) %>% str_replace(pattern=".fastq", replacement=".trimmed.fastq")

    if(ktrim=="l"){ktrim=paste0("ktrim=l")
    } else if (ktrim=="r"){ktrim=paste0("ktrim=r")}

    if(is.numeric(kmer)){kmer=paste0("k=", kmer)
    }  else (kmer=paste0("k=",nchar(sort(primers,decreasing=TRUE)[1])))

    if(is.numeric(maxlength)){maxlength=paste0("maxlength=", maxlength)
    }  else (maxlength="")

    if(is.numeric(mink)){mink=paste0("mink=", mink) # Note - mink makes it noticibly slower
    }  else if (mink==TRUE){mink=paste0("mink=",(ceiling(nchar(sort(primers,decreasing=TRUE)[1])/2)))
    } else if(mink==FALSE){mink=""}


    if(is.numeric(restrictleft)){restrictleft=paste0("restrictleft=",restrictleft)
    }  else (restrictleft=paste0("restrictleft=",nchar(sort(primers)[[1]])))

    if(ordered==TRUE){ordered="ordered=T"} else (ordered="")
    if(is.numeric(hdist)){hdist=paste0("hdist=",hdist)}
    if(copyundefined==TRUE){copyundefined="copyundefined"} else (copyundefined="")
    if(overwrite==TRUE){overwrite="overwrite=TRUE"} else (overwrite="")
    if(tpe==TRUE){tpe="tpe"} else (tpe="")

    if(quality==TRUE){quality="bhist=bhist.txt qhist=qhist.txt gchist=gchist.txt aqhist=aqhist.txt lhist=lhist.txt gcbins=auto"} else (quality="")

    args <- paste(" -cp ",install, in1, in2,literal, restrictleft,out, kmer,mink, hdist,ktrim,tpe, copyundefined, quality, maxlength, overwrite,"-da", collapse = " ")

    print(paste0(args," /n"))
    #Run bbduk
    system2("java", args = args)

  }
  if(nsamples >1){
    for (i in 1:nsamples){
      bbduk(install=install, fwd=fwd[i], rev, primers,restrictleft,outpath, ktrim, ordered, kmer, mink,tpe, hdist,copyundefined, quality, overwrite, maxlength)
    }

  } else if(nsamples == 1){
    bbduk(install, fwd, rev, primers,restrictleft,outpath, ktrim, ordered, kmer, mink,tpe, hdist,copyundefined, overwrite, quality, maxlength)
  }
}



# Demultiplex by primers --------------------------------------------------


#' Demultiplex fusion primers using BBmap Seal
#'
#' @param install
#' @param fwds
#' @param revs
#' @param Fbarcodes
#' @param Rbarcodes
#' @param restrictleft
#' @param outpath
#' @param kmer
#' @param hdist
#' @param copyundefined
#' @param overwrite
#' @param interleaved
#'
#' @return
#' @export
#'
#' @examples
bbtools_demux <- function(install=NULL, fwds, revs=NULL, Fbarcodes=NULL,Rbarcodes,
                          restrictleft=NULL,outpath="trimmed", kmer=NULL, hdist=0,copyundefined=TRUE,
                          overwrite=TRUE, interleaved=FALSE){

  nsamples <- length(fwds)

  bbtools_seal <- function(install=NULL, fwd, rev=NULL, Fbarcodes=NULL,Rbarcodes,
                           restrictleft=NULL,outpath="trimmed", kmer=NULL, hdist=0,copyundefined=TRUE,
                           overwrite=TRUE){

    in1=paste0( "in=",fwd)
    if(!is.null(rev)){in2=paste0( "in2=",rev)}  else (in2="")

    if(!is.null(Fbarcodes) & is.null(Rbarcodes)){
      writeLines(paste0(">Rep", seq(1:length(Fbarcodes))," \n", Fbarcodes),con="Fprimers.fa")
      ref="ref=Fprimers.fa"

    }  else if(!is.null(Fbarcodes) & !is.null(Rbarcodes)){
      writeLines(paste0(">Rep", seq(1:length(Fbarcodes)),"\n", Fbarcodes),con="Fprimers.fa",sep="")
      writeLines(paste0(">Rep", seq(1:length(Rbarcodes)),"\n", Rbarcodes),con="Rprimers.fa",sep="")
      ref="ref=Fprimers.fa,Rprimers.fa"
    }

    pattern = paste0("pattern=",outpath,"/", basename(fwd) %>%
                       str_split_fixed("\\.",n=2) %>%
                       as_tibble() %>% pull(V1) %>%
                       str_replace(pattern="_R1_",replacement="_R1R2_") ,"_%.fastq.gz")

    if(is.numeric(kmer)){kmer=paste0("k=", kmer)
    }  else (kmer=paste0("k=",nchar(sort(c(Fbarcodes,Rbarcodes),decreasing=TRUE)[1])))

    if(is.numeric(restrictleft)){restrictleft=paste0("restrictleft=",restrictleft)
    }  else (restrictleft=paste0("restrictleft=",nchar(sort(c(Fbarcodes,Rbarcodes))[1])))


    if(is.numeric(hdist)){hdist=paste0("hdist=",hdist)
    }
    if(copyundefined==TRUE){copyundefined="copyundefined"} else (copyundefined="")
    if(overwrite==TRUE){overwrite="overwrite=TRUE"} else (overwrite="")

    args <- paste(" -cp ",paste0(install,"/current jgi.Seal"), in1, in2, ref,restrictleft,pattern, kmer,hdist, copyundefined, overwrite, "kpt=t",collapse = " ")

    print(paste0(args," /n"))
    #Run Seal
    system2("java", args = args, stdout = "trim.txt", stderr = "trimerr.txt")

  }


  if(nsamples >1){
    for (i in 1:nsamples){
      bbtools_seal(install=install, fwd=fwds[i], rev=revs[i], Fbarcodes=Fbarcodes,Rbarcodes=Rbarcodes,
                   restrictleft=restrictleft,outpath=outpath, kmer=kmer, hdist=hdist,copyundefined=copyundefined,
                   overwrite=overwrite)
    }

  } else if(nsamples == 1){
    bbtools_seal(install, fwd, rev,
                 Fbarcodes, Rbarcodes,restrictleft,
                 outpath,kmer,hdist,copyundefined,overwrite)
  }
}



# Split interleaved reads -------------------------------------------------


#' Split interleaved reads
#'
#' @param install
#' @param files
#' @param overwrite
#'
#' @return
#' @export
#'
#' @examples
bbtools_split <- function(install=NULL, files, overwrite=FALSE){
  nsamples <- length(files)

  bbtools_reformatreads <- function(install=NULL, file, overwrite=FALSE){
    #Split interleaved reads
    out1 <- paste0("out1=", file %>% str_replace(pattern="_R1R2_",replacement="_R1_"))
    out2 <- paste0("out2=", file %>% str_replace(pattern="_R1R2_",replacement="_R2_"))

    if(overwrite==TRUE){overwrite="overwrite=TRUE"} else (overwrite="")

    reformat_args <- paste(" -cp ",paste0(install,"/current jgi.ReformatReads "),file, out1, out2,overwrite, collapse = " ")
    system2("java", args = reformat_args)
  }

  if(nsamples >1){
    for (i in 1:nsamples){
      bbtools_reformatreads(install=install,file=files[i], overwrite=overwrite)
      file.remove(files[i])
    }

  } else if(nsamples == 1){
    bbtools_reformatreads(install=install,file=files, overwrite=overwrite)
    file.remove(files)
  }
}


