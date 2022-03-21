# Notes ####
# This script is for visulisation of Allosteric Signalling Map (ASM) analysis.
#
# This script is ONLY suitable for SINGLE residue perturbation, and
# to see effects on each functional residue of soluble epoxide hydrolase (sEH).
#
# **For any other uses, please modify accordingly.**
# 
# ************************
# Required inputs:
#   json_all_singular.txt file generated from Python script asm_json.py
#   
# Output:
#   ASM_plot.pdf
#
# ************************

# library ####
{
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(ggplot2)
  library(ggpubr)
  library(ggh4x)
  library(gridExtra) #For grid.arrange
  
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  
}

setwd('../outputs/ASM')

# Input parameters ####
# enter perturbated residue (perturbate one at a time)
selected_resi = c('380','381','382',
                  '507','490','360','359',
                  '523','334','337','338','310','465')

# sequence ####
{
  sequence = 'AACNPSDMSHGYVTVKPRVRLHFVELGSGPAVCLCHGFPESWYSWRYQIPALAQAGYRVLAMDMKGYGESSAPPEIEEYCMEVLCKEMVTFLDKLGLSQAVFIGHDWGGMLVWYMALFYPERVRAVASLNTPFIPANPNMSPLESIKANPVFDYQLYFQEPGVAEAELEQNLSRTFKSLFRASDESVLSMHKVCEAGGLFVNSPEEPSLSRMVTEEEIQFYVQQFKKSGFRGPLNWYRNMERNWKWACKSLGRKILIPALMVTAEKDFVLVPQMSQHMEDWIPHLKRGHIEDCGHWTQMDKPTEVNQILIKWLDSD'
  ol_seq = strsplit(sequence,'')[[1]]
  oneletter = c()
  start = 230
  for (i in (1:length(ol_seq))){
    temp = paste0(ol_seq[i],start)
    start = start +1
    oneletter = append(oneletter,temp)
  }
  oneletter_df = data.frame('resi' = as.character(c(230:545)),'ol_resi' = oneletter)
}


# Read data ####
raw_json_txt = read.table('json_all_singular.txt',header = T, sep = ';',
                            row.names = NULL) %>%
  mutate(`mutation` = as.character(`mutation`),
         `mutation` = gsub('[','',`mutation`,fixed = T),
         `mutation` = gsub(']','',`mutation`,fixed = T)) %>%
  mutate(`domain` = case_when(grepl('383|466',`response`) ~ 'Epoxide Positioners',
                              grepl('335|496|524',`response`) ~'Catalytic Triad',
                              grepl('339|499|336|469',`response`)~'W336 Niche',
                              grepl('417|525|267|408|419|498',`response`) ~ 'F267 Pocket')) %>%
  left_join(oneletter_df,by = c('mutation' = 'resi')) %>%
  rename('old_mutation' = 'mutation') %>%
  rename('mutation' = 'ol_resi') %>%
  mutate(`response` = as.character(`response`)) %>%
  left_join(oneletter_df,by = c('response' = 'resi')) %>%
  select(-response) %>%
  rename('response' = 'ol_resi') %>%
  mutate(`position` = paste0(`domain`,
                             ' (',`response`,')')) 

raw_json_txt$mutation = factor(raw_json_txt$mutation,
                                 levels = oneletter)

raw_json_txt$position = factor(raw_json_txt$position,
                           levels = c('Epoxide Positioners (Y383)','Epoxide Positioners (Y466)',
                                          'Catalytic Triad (D335)','Catalytic Triad (D496)','Catalytic Triad (H524)',
                                          'W336 Niche (W336)','W336 Niche (M339)','W336 Niche (M469)','W336 Niche (L499)',
                                          'F267 Pocket (F267)','F267 Pocket (L408)','F267 Pocket (L417)',
                                          'F267 Pocket (M419)','F267 Pocket (V498)','F267 Pocket (W525)'))



#function ####
gi_plot_func = function(df1,df2,bg_c){
  # df1 = json_plot_A %>%
  #   filter(`position` == 'W336 Niche (M469)')
  if ((all(df1$value)<0)|(all(df1$value)>0)){
    value_list = c(0,df1$value)
  }else{
    value_list = df1$value
  }
  
  y_break_sep = round((max(value_list)-min(value_list))/3,2)
  # y_break = round(c(median(value_list)+y_break_sep,median(value_list),
  #                   median(value_list)-y_break_sep),2)
  y_break = seq(max(value_list),min(value_list),
                by = (min(value_list)-max(value_list))/3)
  y_label = format(round(y_break,1),nsmall = 1)
  
  fig =   ggplot(data = df1, aes(x = mutation, y = value))+
    geom_point(fill = 'grey90',shape = 21,size = 7) +
    geom_hline(yintercept = 0, linetype = 'dashed', col = 'red',size = 1)+
    geom_point(data = df2,aes(x = mutation,y = value),fill = 'red',
               shape = 21, size = 7)+
    scale_y_continuous(breaks = y_break,labels = y_label)+
    theme(axis.text.x = element_text(size = 24,face = 'bold'),
      axis.text.y = element_text(size = 24, face = 'bold'),
      strip.text = element_text(size = 24,face = 'bold',
                                margin = margin(0.3,0,0.3,0, "cm")),
      strip.background = element_rect(fill = bg_c,color = 'black',size = 1),
      panel.grid = element_line(color = 'grey20'),
      panel.grid.major = element_line(color = 'grey20'),
      panel.background = element_rect(color = 'black',size = 1,fill = 'white'))+
    labs(x = '',
         y = '')+
    facet_wrap(~position,ncol = 1)
  return(fig)
}

# gi plot ####
all_resi = c('Epoxide Positioners (Y383)','Epoxide Positioners (Y466)',
             'Catalytic Triad (D335)','Catalytic Triad (D496)','Catalytic Triad (H524)',
             'W336 Niche (W336)','W336 Niche (M339)','W336 Niche (M469)','W336 Niche (L499)',
             'F267 Pocket (F267)','F267 Pocket (L408)','F267 Pocket (L417)',
             'F267 Pocket (M419)','F267 Pocket (V498)','F267 Pocket (W525)')

for (i in (1:length(all_resi))){
  r = all_resi[i]
  if (grepl('Epoxide Positioners',r)){
    cl = '#9ac5f0' #marine
    name = paste('ep',i,sep = '_')
  } else{
    if (grepl('Catalytic Triad',r)){
      cl = '#a1a1d8' #bluepurple
      name = paste('ct',i,sep = '_')
    } else{
      if (grepl('W336 Niche',r)){
        cl = '#73b9b9' #cyan
        name = paste('w3',i,sep = '_')
      }else{
        cl = '#dfa0df'
        name = paste('f2',i,sep = '_')
      }
    }
  }
  plot_all = raw_json_txt %>%
    filter(r == `position`) %>%
    filter(old_mutation %in% selected_resi)
  
  plot_hl = plot_all %>%
    filter(abs(`value`)>=2)
  
  current_plot = gi_plot_func(plot_all,plot_hl,cl)
  assign(name,current_plot)
  
}

all_p = grid.arrange(ep_1,ep_2,ct_3,ct_4,ct_5,
                     w3_6,w3_7,w3_8,w3_9,
                     f2_10,f2_11,f2_12,f2_13,f2_14,f2_15,ncol = 3)

assign(paste('plot','others',sep = '_'),all_p)


#Export ####

ggpubr::ggarrange(
  plot_others,
  nrow = 1, ncol = 1
) %>% #IMPORTANT to have this chain ##
  ggexport(filename = "ASM_plot.pdf", width = 42, height = 12)

