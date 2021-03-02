require(rmarkdown)

rmarkdown::render('R/Preregistration.Rmd', output_dir = 'Protocols')
#rmarkdown::render('R/Preregistration.Rmd', output_dir = 'Protocols', word_document(toc = TRUE))
#rmarkdown::render('R/Preregistration.Rmd', output_dir = 'Protocols', word_document(highlight = "zenburn"))

rmarkdown::render('R/protocol_sperm.Rmd', output_dir = 'Protocols')
#rmarkdown::render('R/protocol_sperm.Rmd', output_dir = 'Protocols', word_document(toc = TRUE))