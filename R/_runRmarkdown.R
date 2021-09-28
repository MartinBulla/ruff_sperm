require(rmarkdown)

rmarkdown::render('R/motility.R', output_dir = 'Output')
rmarkdown::render('R/EXP_test-40.Rmd', output_dir = 'Output')

rmarkdown::render('R/Methods.Rmd', output_dir = 'Output')

rmarkdown::render('R/Preregistration_v3.Rmd', output_dir = 'Protocols')

rmarkdown::render('R/Preregistration_v2.Rmd', output_dir = 'Protocols')

rmarkdown::render('R/Preregistration.Rmd', output_dir = 'Protocols')
#rmarkdown::render('R/Preregistration.Rmd', output_dir = 'Protocols', word_document(toc = TRUE))
#rmarkdown::render('R/Preregistration.Rmd', output_dir = 'Protocols', word_document(highlight = "zenburn"))

rmarkdown::render('R/protocol_sperm.Rmd', output_dir = 'Protocols')
#rmarkdown::render('R/protocol_sperm.Rmd', output_dir = 'Protocols', word_document(toc = TRUE))