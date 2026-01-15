getwd()

# Antes eu tive que executar esse comando no terminal para tirar a vigilancia do git da minha pasta princiapl rm -rf ~/.git

library(usethis)

use_git_config(user.name = "Seu Nome", user.email = "seu-email@exemplo.com")

gitcreds::gitcreds_set() # ver token

usethis::use_git() # responda sim para commitment

usethis::use_github() # cria automaticamente o github no repositorio online

# Depois vou em version control para habilitar o R


