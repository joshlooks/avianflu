# Define a transition function for each transition
# T->E->I-> (upper)
# T->E->I-> (lower)
# Production, clearance, infection, transport, transport
# Production, clearance, infection, immune
# Adaptive, growth
transitions_trans = list(c(Tu = -1, Eu = +1), #Tu -> Eu
                   c(Eu = -1, Iu = +1), #Eu -> Iu
                   c(Iu = -1, Du = +1),          #Iu -> D
                   c(Tl = -1, El = +1), #Tl -> El
                   c(El = -1, Il = +1), #El -> Il
                   c(Il = -1, Dl = +1),          #Il -> D
                   c(Vu = +1),          #Production from Iu
                   c(Vu = -1),          #Clearance of Vu
                   c(Vu = -1),          #Infection from Vu
                   c(Vu = -1, Vl = +1), #Vu -> Vl transport
                   c(Vu = +1, Vl = -1), #Vl -> Vu transport
                   c(Vl = +1),          #Production from Il
                   c(Vl = -1),          #Clearance of Vl
                   c(Vl = -1),          #Infection from Vl
                   c(Vl = -1),          #Immune response to Vl
                   c(X = +1),           #Adaptive growth from Vl
                   c(X = +1))           #Adaptive exponential growth
# Define a rate function for:
# T->E->I (upper)
# T->E->I (lower)
# Production, clearance, infection, transport, transport
# Production, clearance, infection, immune
# Adaptive, growth
lvrates_trans <- function(x, para, t){
  return(c(infu=para$b_u*x['Tu']*x['Vu'],             #Tu -> Eu
           onu=para$g*x['Eu'],                        #Eu -> Iu
           Du=para$d*x['Iu'],                         #Iu -> D
           infl=para$b_l*x['Tl']*x['Vl'],             #Tl -> El
           onl=para$g*x['El'],                        #El -> Il
           Dl=para$d*x['Il'],                         #Il -> D
           pru=para$pu*x['Iu'],                       #Production from Iu
           clu=para$c*x['Vu'],                        #Clearance of Vu
           infvu=para$gamma*para$b_u*x['Tu']*x['Vu'], #Infection from Vu
           tu=para$t_1*x['Vu'],                       #Vu -> Vl transport
           tl=para$t_2*x['Vl'],                       #Vl -> Vu transport
           prl=para$pl*x['Il'],                       #Production from Il
           cll=para$c*x['Vl'],                        #Clearance of Vl
           infvl=para$gamma*para$b_l*x['Tl']*x['Vl'], #Infection from Vl
           iml=para$k*x['Vl']*x['X'],                 #Immune response to Vl
           adx=para$f*x['Vl'],                        #Adaptive growth from Vl
           gwx=para$r*x['X']))                        #Adaptive exponential growth
}

# Define a transition function for each transition
# T->E->I-> (upper)
# T->E->I-> (lower)
# Production, clearance, infection, diffusion, diffusion, advection
# Production, clearance, infection, immune
# Adaptive, growth
transitions_diffadv = list(c(Tu = -1, Eu = +1), #Tu -> Eu
                   c(Eu = -1, Iu = +1), #Eu -> Iu
                   c(Iu = -1, Du = +1),          #Iu -> D
                   c(Tl = -1, El = +1), #Tl -> El
                   c(El = -1, Il = +1), #El -> Il
                   c(Il = -1, Dl = +1),          #Il -> D
                   c(Vu = +1),          #Production from Iu
                   c(Vu = -1),          #Clearance of Vu
                   c(Vu = -1),          #Infection from Vu
                   c(Vu = -1, Vl = +1), #Vu -> Vl diffusion
                   c(Vu = +1, Vl = -1), #Vl -> Vu diffusion
                   c(Vu = +1, Vl = -1), #Vl -> Vu advection
                   c(Vl = +1),          #Production from Il
                   c(Vl = -1),          #Clearance of Vl
                   c(Vl = -1),          #Infection from Vl
                   c(Vl = -1),          #Immune response to Vl
                   c(X = +1),           #Adaptive growth from Vl
                   c(X = +1))           #Adaptive exponential growth
# Define a rate function for:
# T->E->I (upper)
# T->E->I (lower)
# Production, clearance, infection, diffusion, diffusion, advection
# Production, clearance, infection, immune
# Adaptive, growth
lvrates_diffadv <- function(x, para, t){
  return(c(infu=para$b_u*x['Tu']*x['Vu'],             #Tu -> Eu
           onu=para$g*x['Eu'],                        #Eu -> Iu
           Du=para$d*x['Iu'],                         #Iu -> D
           infl=para$b_l*x['Tl']*x['Vl'],             #Tl -> El
           onl=para$g*x['El'],                        #El -> Il
           Dl=para$d*x['Il'],                         #Il -> D
           pru=para$pu*x['Iu'],                       #Production from Iu
           clu=para$c*x['Vu'],                        #Clearance of Vu
           infvu=para$gamma*para$b_u*x['Tu']*x['Vu'], #Infection from Vu
           Diffu=max(0, para$D*(x['Vu']-x['Vl'])),    #Vu -> Vl diffusion
           Diffl=max(0, para$D*(x['Vl']-x['Vu'])),    #Vl -> Vu diffusion
           Adv=para$a*x['Vl'],                        #Vl -> Vu advection
           prl=para$pl*x['Il'],                       #Production from Il
           cll=para$c*x['Vl'],                        #Clearance of Vl
           infvl=para$gamma*para$b_l*x['Tl']*x['Vl'], #Infection from Vl
           iml=para$k*x['Vl']*x['X'],                 #Immune response to Vl
           adx=para$f*x['Vl'],                        #Adaptive growth from Vl
           gwx=para$r*x['X']))                        #Adaptive exponential growth
}
