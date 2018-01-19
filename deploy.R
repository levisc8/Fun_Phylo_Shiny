# Deploy Script

# NOTE: Before deploying, check versions of mgcv on machine and at
# http://docs.rstudio.com/shinyapps.io/appendix.html#default-system-packages.
# If local version is more recent, downgrade to version ^^ and run deployment


rsconnect::deployApp(appFiles = c('ui.R' , 'server.R', 'Readme_general.md'),
                     logLevel = 'verbose')
