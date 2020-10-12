# import library
library(shiny)
library(shinythemes)
library(tidyverse)
# grid.arrange function
library(gridExtra)

# import formula
source('formula.R')

ui<-navbarPage(title='Model Type',
               theme=shinytheme('united'),
               tabPanel('Black Scholes Model',
                        headerPanel('Black Scholes Model, Binomial Tree Model, and Monte Carlo Simulation'),
                        sidebarLayout(
                            sidebarPanel(width=3,
                                         textInput(inputId='tab1_S',
                                                   label='Underlying Price (S):',
                                                   value=100),
                                         textInput(inputId='tab1_K',
                                                   label='Strike Price (K):',
                                                   value=70),
                                         textInput(inputId='tab1_r',
                                                   label='Risk Free Interest Rate (r):',
                                                   value=0.01),
                                         textInput(inputId='tab1_maturity',
                                                   label='Time to Maturity (T):',
                                                   value=1),
                                         textInput(inputId='tab1_sig',
                                                   label=HTML('Volatility of Underlying Return (&sigma;):'),
                                                   value=0.2),
                                         textInput(inputId='tab1_tree_n',
                                                   label='Number of Period for Binomial Tree Model :',
                                                   value=100),
                                         textInput(inputId='tab1_sim_n',
                                                   label='Number of Path for Monte Carlo Simulation :',
                                                   value=10000),
                                         br(),
                                         h4('Ignore Dividend Payment')),
                            mainPanel(tabsetPanel(type='tabs',
                                                  tabPanel('Price',
                                                           h3('Price'),
                                                           br(),
                                                           h4('Black Scholes Model'),
                                                           tableOutput('bs_price_table'),
                                                           h4('Binomial Tree Model (European)'),
                                                           tableOutput('bs_bt_euro_price_table'),
                                                           h4('Binomial Tree Model (American)'),
                                                           tableOutput('bs_bt_amer_price_table'),
                                                           h4('Monte Carlo Simulation'),
                                                           tableOutput('bs_mc_price_table'),
                                                           h3('Other Information'),
                                                           br(),
                                                           h4('Black Scholes Model'),
                                                           tableOutput('bs_N_d_table'),
                                                           tableOutput('bs_greek_table'),
                                                           h4('Binomial Tree Model'),
                                                           tableOutput('bs_bt_arg_table'),
                                                           h4('Monte Carlo Simulation'),
                                                           tableOutput('bs_mc_ste_table')),
                                                  tabPanel('Plot',
                                                           h3('Convergence of European Option Price'),
                                                           plotOutput('bs_bt_euro_cov_plot',width='100%'),
                                                           plotOutput('bs_mc_euro_cov_plot',width='100%')))))),
                 tabPanel('Black Scholes Merton Model',
                          headerPanel('Black Scholes Merton Model, Binomial Tree Model and Monte Carlo Simulation'),
                          sidebarLayout(
                            sidebarPanel(width=3,
                                         textInput(inputId='tab2_S',
                                                   label='Underlying Price (S):',
                                                   value=100),
                                         textInput(inputId='tab2_K',
                                                   label='Strike Price (K):',
                                                   value=70),
                                         textInput(inputId='tab2_r',
                                                   label='Risk Free Interest Rate (r):',
                                                   value=0.01),
                                         textInput(inputId='tab2_q',
                                                   label='Dividend Yield (q):',
                                                   value=0.02),
                                         textInput(inputId='tab2_maturity',
                                                   label='Time to Maturity (T):',
                                                   value=1),
                                         textInput(inputId='tab2_sig',
                                                   label=HTML('Volatility of Underlying Return (&sigma;):'),
                                                   value=0.2),
                                         textInput(inputId='tab2_tree_n',
                                                   label='Number of Period for Binomial Tree Model :',
                                                   value=100),
                                         textInput(inputId='tab2_sim_n',
                                                   label='Number of Path for Monte Carlo Simulation :',
                                                   value=10000),
                                         br(),
                                         h4('With Dividend Payment')),
                            mainPanel(tabsetPanel(type='tabs',
                                                  tabPanel('Price',
                                                           h3('Price'),
                                                           br(),
                                                           h4('Black Scholes Merton Model'),
                                                           tableOutput('bsm_price_table'),
                                                           h4('Binomial Tree Model (European)'),
                                                           tableOutput('bsm_bt_euro_price_table'),
                                                           h4('Binomial Tree Model (American)'),
                                                           tableOutput('bsm_bt_amer_price_table'),
                                                           h4('Monte Carlo Simulation'),
                                                           tableOutput('bsm_mc_price_table'),
                                                           h3('Other Information'),
                                                           br(),
                                                           h4('Black Scholes Merton Model'),
                                                           tableOutput('bsm_N_d_table'),
                                                           tableOutput('bsm_greek_table'),
                                                           h4('Binomial Tree Model'),
                                                           tableOutput('bsm_bt_arg_table'),
                                                           h4('Monte Carlo Simulation'),
                                                           tableOutput('bsm_mc_ste_table')),
                                                  tabPanel('Plot',
                                                           h3('Convergence of European Option Price'),
                                                           plotOutput('bsm_bt_euro_cov_plot',width='100%'),
                                                           plotOutput('bsm_mc_euro_cov_plot',width='100%'))))
                            )),
                 tabPanel('Black Model',
                          headerPanel('Black Model, Binomial Tree Model and Monte Carlo Simulation'),
                          sidebarLayout(
                            sidebarPanel(width=3,
                                         textInput(inputId='tab3_F',
                                                   label='Forward Price (F):',
                                                   value=100),
                                         textInput(inputId='tab3_K',
                                                   label='Strike Price (K):',
                                                   value=70),
                                         textInput(inputId='tab3_r',
                                                   label='Risk Free Interest Rate (r):',
                                                   value=0.01),
                                         textInput(inputId='tab3_maturity',
                                                   label='Time to Maturity (T):',
                                                   value=1),
                                         textInput(inputId='tab3_sig',
                                                   label=HTML('Volatility of Underlying Return (&sigma;):'),
                                                   value=0.2),
                                         textInput(inputId='tab3_tree_n',
                                                   label='Number of Period for Binomial Tree Model :',
                                                   value=100),
                                         textInput(inputId='tab3_sim_n',
                                                   label='Number of Path for Monte Carlo Simulation :',
                                                   value=10000),
                                         br()),
                            mainPanel(tabsetPanel(type='tabs',
                                                  tabPanel('Price',
                                                           h3('Price'),
                                                           br(),
                                                           h4('Black Model'),
                                                           tableOutput('bl_price_table'),
                                                           h4('Binomial Tree Model (European)'),
                                                           tableOutput('bl_bt_euro_price_table'),
                                                           h4('Binomial Tree Model (American)'),
                                                           tableOutput('bl_bt_amer_price_table'),
                                                           h4('Monte Carlo Simulation'),
                                                           tableOutput('bl_mc_price_table'),
                                                           h3('Other Information'),
                                                           br(),
                                                           h4('Black Model'),
                                                           tableOutput('bl_N_d_table'),
                                                           tableOutput('bl_greek_table'),
                                                           h4('Binomial Tree Model'),
                                                           tableOutput('bl_bt_arg_table'),
                                                           h4('Monte Carlo Simulation'),
                                                           tableOutput('bl_mc_ste_table')),
                                                  tabPanel('Plot',
                                                           h3('Convergence of European Option Price'),
                                                           plotOutput('bl_bt_euro_cov_plot',width='100%'),
                                                           plotOutput('bl_mc_euro_cov_plot',width='100%'))))
                            ))
               )

# Define server logic required to draw a histogram
server<-function(input,output){
  # tab1
  # Dataframe
  bs_price_df<-reactive({
    S<-as.numeric(input$tab1_S)
    K<-as.numeric(input$tab1_K)
    r<-as.numeric(input$tab1_r)
    maturity<-as.numeric(input$tab1_maturity)
    sig<-as.numeric(input$tab1_sig)
    
    make_bs_price_df(S,K,r,maturity,sig)
  })
  
  
  bs_bt_euro_price_df<-reactive({
    S<-as.numeric(input$tab1_S)
    K<-as.numeric(input$tab1_K)
    r<-as.numeric(input$tab1_r)
    maturity<-as.numeric(input$tab1_maturity)
    sig<-as.numeric(input$tab1_sig)
    tree_n<-as.numeric(input$tab1_tree_n)
    
    make_bs_bt_euro_price_df(S,K,r,maturity,sig,tree_n)
  })
  
  
  bs_bt_amer_price_df<-reactive({
    S<-as.numeric(input$tab1_S)
    K<-as.numeric(input$tab1_K)
    r<-as.numeric(input$tab1_r)
    maturity<-as.numeric(input$tab1_maturity)
    sig<-as.numeric(input$tab1_sig)
    tree_n<-as.numeric(input$tab1_tree_n)
    
    make_bs_bt_amer_price_df(S,K,r,maturity,sig,tree_n)
  })
  
  
  bs_mc_price_df<-reactive({
    S<-as.numeric(input$tab1_S)
    K<-as.numeric(input$tab1_K)
    r<-as.numeric(input$tab1_r)
    maturity<-as.numeric(input$tab1_maturity)
    sig<-as.numeric(input$tab1_sig)
    sim_n<-as.numeric(input$tab1_sim_n)
    
    make_bs_mc_price_df(S,K,r,maturity,sig,sim_n)
  })
  
  
  # tab2
  # Dataframe
  bsm_price_df<-reactive({
    S<-as.numeric(input$tab2_S)
    K<-as.numeric(input$tab2_K)
    r<-as.numeric(input$tab2_r)
    q<-as.numeric(input$tab2_q)
    maturity<-as.numeric(input$tab2_maturity)
    sig<-as.numeric(input$tab2_sig)
    
    make_bsm_price_df(S,K,r,q,maturity,sig)
  })
  
  
  bsm_bt_euro_price_df<-reactive({
    S<-as.numeric(input$tab2_S)
    K<-as.numeric(input$tab2_K)
    r<-as.numeric(input$tab2_r)
    q<-as.numeric(input$tab2_q)
    maturity<-as.numeric(input$tab2_maturity)
    sig<-as.numeric(input$tab2_sig)
    tree_n<-as.numeric(input$tab2_tree_n)
    
    make_bsm_bt_euro_price_df(S,K,r,q,maturity,sig,tree_n)
  })
  
  
  bsm_bt_amer_price_df<-reactive({
    S<-as.numeric(input$tab2_S)
    K<-as.numeric(input$tab2_K)
    r<-as.numeric(input$tab2_r)
    q<-as.numeric(input$tab2_q)
    maturity<-as.numeric(input$tab2_maturity)
    sig<-as.numeric(input$tab2_sig)
    tree_n<-as.numeric(input$tab2_tree_n)
    
    make_bsm_bt_amer_price_df(S,K,r,q,maturity,sig,tree_n)
  })
  
  
  bsm_mc_price_df<-reactive({
    S<-as.numeric(input$tab2_S)
    K<-as.numeric(input$tab2_K)
    r<-as.numeric(input$tab2_r)
    q<-as.numeric(input$tab2_q)
    maturity<-as.numeric(input$tab2_maturity)
    sig<-as.numeric(input$tab2_sig)
    sim_n<-as.numeric(input$tab2_sim_n)
    
    make_bsm_mc_price_df(S,K,r,q,maturity,sig,sim_n)
  })
  
  
  # tab3
  # Dataframe
  bl_price_df<-reactive({
    Forward<-as.numeric(input$tab3_F)
    K<-as.numeric(input$tab3_K)
    r<-as.numeric(input$tab3_r)
    maturity<-as.numeric(input$tab3_maturity)
    sig<-as.numeric(input$tab3_sig)
    
    make_bl_price_df(Forward,K,r,maturity,sig)
  })
  
  
  bl_bt_euro_price_df<-reactive({
    Forward<-as.numeric(input$tab3_F)
    K<-as.numeric(input$tab3_K)
    r<-as.numeric(input$tab3_r)
    maturity<-as.numeric(input$tab3_maturity)
    sig<-as.numeric(input$tab3_sig)
    tree_n<-as.numeric(input$tab3_tree_n)
    
    make_bl_bt_euro_price_df(Forward,K,r,maturity,sig,tree_n)
  })
  
  
  bl_bt_amer_price_df<-reactive({
    Forward<-as.numeric(input$tab3_F)
    K<-as.numeric(input$tab3_K)
    r<-as.numeric(input$tab3_r)
    maturity<-as.numeric(input$tab3_maturity)
    sig<-as.numeric(input$tab3_sig)
    tree_n<-as.numeric(input$tab3_tree_n)
    
    make_bl_bt_amer_price_df(Forward,K,r,maturity,sig,tree_n)
  })
  
  
  bl_mc_price_df<-reactive({
    Forward<-as.numeric(input$tab3_F)
    K<-as.numeric(input$tab3_K)
    r<-as.numeric(input$tab3_r)
    maturity<-as.numeric(input$tab3_maturity)
    sig<-as.numeric(input$tab3_sig)
    sim_n<-as.numeric(input$tab3_sim_n)
    
    make_bl_mc_price_df(Forward,K,r,maturity,sig,sim_n)
  })
  
  
  
  # tab1
  # Plot
  output$bs_bt_euro_cov_plot<-renderPlot({
    S<-as.numeric(input$tab1_S)
    K<-as.numeric(input$tab1_K)
    r<-as.numeric(input$tab1_r)
    maturity<-as.numeric(input$tab1_maturity)
    sig<-as.numeric(input$tab1_sig)
    tree_n<-as.numeric(input$tab1_tree_n)
    
    plot_bs_bt_euro_cov(S,K,r,maturity,sig,tree_n)
  })
  
  
  output$bs_mc_euro_cov_plot<-renderPlot({
    S<-as.numeric(input$tab1_S)
    K<-as.numeric(input$tab1_K)
    r<-as.numeric(input$tab1_r)
    maturity<-as.numeric(input$tab1_maturity)
    sig<-as.numeric(input$tab1_sig)
    sim_n<-as.numeric(input$tab1_sim_n)
    
    plot_bs_mc_euro_cov(S,K,r,maturity,sig,sim_n)
  })
  
  
  
  # tab2
  # Plot
  output$bsm_bt_euro_cov_plot<-renderPlot({
    S<-as.numeric(input$tab2_S)
    K<-as.numeric(input$tab2_K)
    r<-as.numeric(input$tab2_r)
    q<-as.numeric(input$tab2_q)
    maturity<-as.numeric(input$tab2_maturity)
    sig<-as.numeric(input$tab2_sig)
    tree_n<-as.numeric(input$tab2_tree_n)
    
    plot_bsm_bt_euro_cov(S,K,r,q,maturity,sig,tree_n)
  })
  
  
  output$bsm_mc_euro_cov_plot<-renderPlot({
    S<-as.numeric(input$tab2_S)
    K<-as.numeric(input$tab2_K)
    r<-as.numeric(input$tab2_r)
    q<-as.numeric(input$tab2_q)
    maturity<-as.numeric(input$tab2_maturity)
    sig<-as.numeric(input$tab2_sig)
    sim_n<-as.numeric(input$tab2_sim_n)
    
    plot_bsm_mc_euro_cov(S,K,r,q,maturity,sig,sim_n)
  })
  
  
  
  # tab3
  # Plot
  output$bl_bt_euro_cov_plot<-renderPlot({
    Forward<-as.numeric(input$tab3_F)
    K<-as.numeric(input$tab3_K)
    r<-as.numeric(input$tab3_r)
    maturity<-as.numeric(input$tab3_maturity)
    sig<-as.numeric(input$tab3_sig)
    tree_n<-as.numeric(input$tab3_tree_n)
    
    plot_bl_bt_euro_cov(Forward,K,r,maturity,sig,tree_n)
  })
  
  
  output$bl_mc_euro_cov_plot<-renderPlot({
    Forward<-as.numeric(input$tab3_F)
    K<-as.numeric(input$tab3_K)
    r<-as.numeric(input$tab3_r)
    maturity<-as.numeric(input$tab3_maturity)
    sig<-as.numeric(input$tab3_sig)
    sim_n<-as.numeric(input$tab3_sim_n)
    
    plot_bl_mc_euro_cov(Forward,K,r,maturity,sig,sim_n)
  })

  
  
  # tab1
  output$bs_price_table<-renderTable({
    bs_price_df()[[1]]},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  output$bs_N_d_table<-renderTable({
    bs_price_df()[[2]]},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  output$bs_greek_table<-renderTable({
    bs_price_df()[[3]]},
    rownames=TRUE,digits=4,width='50%',align='c')
  
  output$bs_bt_euro_price_table<-renderTable({
    bs_bt_euro_price_df()},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  output$bs_bt_amer_price_table<-renderTable({
    bs_bt_amer_price_df()[[1]]},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  output$bs_bt_arg_table<-renderTable({
    bs_bt_amer_price_df()[[2]]},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  output$bs_mc_price_table<-renderTable({
    bs_mc_price_df()[[1]]},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  output$bs_mc_ste_table<-renderTable({
    bs_mc_price_df()[[2]]},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  
  # tab2
  output$bsm_price_table<-renderTable({
    bsm_price_df()[[1]]},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  output$bsm_N_d_table<-renderTable({
    bsm_price_df()[[2]]},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  output$bsm_greek_table<-renderTable({
    bsm_price_df()[[3]]},
    rownames=TRUE,digits=4,width='50%',align='c')
  
  output$bsm_bt_euro_price_table<-renderTable({
    bsm_bt_euro_price_df()},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  output$bsm_bt_amer_price_table<-renderTable({
    bsm_bt_amer_price_df()[[1]]},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  output$bsm_bt_arg_table<-renderTable({
    bsm_bt_amer_price_df()[[2]]},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  output$bsm_mc_price_table<-renderTable({
    bsm_mc_price_df()[[1]]},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  output$bsm_mc_ste_table<-renderTable({
    bsm_mc_price_df()[[2]]},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  
  # tab3
  output$bl_price_table<-renderTable({
    bl_price_df()[[1]]},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  output$bl_N_d_table<-renderTable({
    bl_price_df()[[2]]},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  output$bl_greek_table<-renderTable({
    bl_price_df()[[3]]},
    rownames=TRUE,digits=4,width='50%',align='c')
  
  output$bl_bt_euro_price_table<-renderTable({
    bl_bt_euro_price_df()},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  output$bl_bt_amer_price_table<-renderTable({
    bl_bt_amer_price_df()[[1]]},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  output$bl_bt_arg_table<-renderTable({
    bl_bt_amer_price_df()[[2]]},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  output$bl_mc_price_table<-renderTable({
    bl_mc_price_df()[[1]]},
    rownames=FALSE,digits=4,width='50%',align='c')
  
  output$bl_mc_ste_table<-renderTable({
    bl_mc_price_df()[[2]]},
    rownames=FALSE,digits=4,width='50%',align='c')
}

# Run the application 
shinyApp(ui=ui,server=server)
