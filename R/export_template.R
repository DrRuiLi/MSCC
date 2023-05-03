

export_template <- function( export.dir = choose.files() ){

  write.xlsx(compound.template,file = export.dir)


}
