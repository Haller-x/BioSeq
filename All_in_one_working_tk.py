####DNA COMPLEMENTAR###
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SubsMat import MatrixInfo as matlist
from tkinter import filedialog
import tkinter.messagebox
import tkinter as tk
import pyperclip
from tkinter import ttk
from dict_values import *
'''
alterações a serem feitas:
ADICIONAR/ARRUMAR A FUNÇÃO DNA SEQ QUE PREENCHE AS ENTRIES (E BOTAO)
ADICIONAR/ARRUMAR A FUNCAO PAIRWISE-NT E PROT ( E BOTAO)
PERMITIR SALVAR ALINHAMENTOS EM ARQUIVO E/OU GERAR AVISO COM O SCORE
REALIZAR O TRATAMENTO DE EXECEÇÃO
'''
root= tk.Tk()
root.title('BioSeq Tool - by Haller-x (GLA7)')
canvas1 = tk.Canvas(root, width = 600, height = 600,  relief = 'raised')
canvas1.pack()

#labels
label_bio_seq = tk.Label(root, text='Bio Seq')
label_bio_seq.config(font=('helvetica', 25,'bold'))
canvas1.create_window(300, 50, window=label_bio_seq)

label_sub = tk.Label(root, text='Transcription, Translation and Replication')
label_sub.config(font=('helvetica', 12,'bold'))
canvas1.create_window(300, 100, window=label_sub)

####labels Transcription, Translation and Replication
label_dna_seq = tk.Label(root, text='DNA seq:')
label_dna_seq.config(font=('helvetica', 11))
canvas1.create_window(195, 140, window=label_dna_seq)#


label_complementary_dna = tk.Label(root, text='Complementary DNA:')
label_complementary_dna.config(font=('helvetica', 11))
canvas1.create_window(155, 180, window=label_complementary_dna)


label_rna = tk.Label(root, text='RNA:')
label_rna.config(font=('helvetica', 11))
canvas1.create_window(210, 220, window=label_rna)


label_seq_aa = tk.Label(root, text='Seq AA:')
label_seq_aa.config(font=('helvetica', 11))
canvas1.create_window(200, 260, window=label_seq_aa)

###need to change it
label__ = tk.Label(root, text='_______________________________________________________________________________________________________________')###
label__.config(font=('helvetica', 12))
canvas1.create_window(100, 300, window=label__)

#labels pairwise
label_pairwise = tk.Label(root, text='  Pairwise simple alignment-Nucleotides')###
label_pairwise.config(font=('helvetica', 11,'bold'))
canvas1.create_window(140, 340, window=label_pairwise)

label_dna_seq1_pairwise = tk.Label(root, text='First Seq:')###
label_dna_seq1_pairwise.config(font=('helvetica', 11))
canvas1.create_window(80, 370, window=label_dna_seq1_pairwise)

label_dna_seq2_pairwise = tk.Label(root, text='Second Seq:')###
label_dna_seq2_pairwise.config(font=('helvetica', 11))
canvas1.create_window(70, 400, window=label_dna_seq2_pairwise)

label_match_weight_pairwise = tk.Label(root, text='Match weight:')###
label_match_weight_pairwise.config(font=('helvetica', 11))
canvas1.create_window(70, 430, window=label_match_weight_pairwise)

label_mismatch_weight_pairwise = tk.Label(root, text='Mismatch weight:')###
label_mismatch_weight_pairwise.config(font=('helvetica', 11))
canvas1.create_window(70, 460, window=label_mismatch_weight_pairwise)

label_gap_open = tk.Label(root, text='Gap open:')###
label_gap_open.config(font=('helvetica', 11))
canvas1.create_window(70, 490, window=label_gap_open)

label_gap_extend = tk.Label(root, text='Gap extend:')###
label_gap_extend.config(font=('helvetica', 11))
canvas1.create_window(70, 520, window=label_gap_extend)

#labels matrix alignment

label_pairwise_prot = tk.Label(root, text='  Pairwise simple alignment-Protein')###
label_pairwise_prot.config(font=('helvetica', 11,'bold'))
canvas1.create_window(425, 340, window=label_pairwise_prot)

label_dna_seq1_matrix = tk.Label(root, text='First protein:')###
label_dna_seq1_matrix.config(font=('helvetica', 11))
canvas1.create_window(380, 370, window=label_dna_seq1_matrix)

label_dna_seq2_matrix = tk.Label(root, text='Second protein:')###
label_dna_seq2_matrix.config(font=('helvetica', 11))
canvas1.create_window(370, 400, window=label_dna_seq2_matrix)

label_matrix = tk.Label(root, text='Matrix:')###
label_matrix.config(font=('helvetica', 11))
canvas1.create_window(395, 430, window=label_matrix)
var = tk.StringVar(root)
var.set("matlist.benner6")
optionMatrix = tk.OptionMenu(root, var, 'matlist.benner6', 'matlist.benner22', 'matlist.benner74', 'matlist.blosum100', 'matlist.blosum30', 'matlist.blosum35', 'matlist.blosum40', 'matlist.blosum45', 'matlist.blosum50','matlist.blosum55', 'matlist.blosum60', 'matlist.blosum62', 'matlist.blosum65', 'matlist.blosum70','matlist.blosum75', 'matlist.blosum80', 'matlist.blosum85', 'matlist.blosum90', 'matlist.blosum95','matlist.feng', 'matlist.fitch', 'matlist.genetic', 'matlist.gonnet', 'matlist.grant', 'matlist.ident', 'matlist.johnson', 'matlist.levin', 'matlist.mclach','matlist.miyata', 'matlist.nwsgappep', 'matlist.pam120', 'matlist.pam180', 'matlist.pam250', 'matlist.pam30', 'matlist.pam300', 'matlist.pam60', 'matlist.pam90', 'matlist.rao', 'matlist.risler', 'matlist.structure')
canvas1.create_window(485, 430, window=optionMatrix)
chosen_matrix = var.get()
chosen = dict_matrix[chosen_matrix]


#entries DNA PART

entry_dna_seq = tk.Entry (root)#dna seq
canvas1.create_window(300, 140, window=entry_dna_seq)

entry_complementary_dna = tk.Entry (root)#complementary dna
canvas1.create_window(300, 180, window=entry_complementary_dna)

entry_rna = tk.Entry (root)#rna
canvas1.create_window(300, 220, window=entry_rna)

entry_seq_aa = tk.Entry (root)#seq aa
canvas1.create_window(300, 260, window=entry_seq_aa)

#entries alignment part
entry_seq_dna1_pairwise = tk.Entry (root)
canvas1.create_window(195, 370, window=entry_seq_dna1_pairwise)

entry_seq_dna2_pairwise = tk.Entry (root)
canvas1.create_window(195, 400, window=entry_seq_dna2_pairwise)


entry_match_weight_pairwise = tk.Entry (root)
canvas1.create_window(195, 430, window=entry_match_weight_pairwise)


entry_mismatch_weight_pairwise = tk.Entry (root)
canvas1.create_window(195, 460, window=entry_mismatch_weight_pairwise)


entry_gap_open = tk.Entry (root)#seq aa
canvas1.create_window(195, 490, window=entry_gap_open)


entry_gap_extend = tk.Entry (root)#seq aa
canvas1.create_window(195, 520, window=entry_gap_extend)


#entries matrix alignment

entry_dna_seq1_matrix = tk.Entry (root)#seq aa
canvas1.create_window(485, 370, window=entry_dna_seq1_matrix)


entry_dna_seq2_matrix = tk.Entry (root)#seq aa
canvas1.create_window(485, 400, window=entry_dna_seq2_matrix)

#functions TTR
#notworking yet

def generatedna(molde):
        dnaM = Seq(molde, generic_dna)

        #Fita complementar
        dnaC = dnaM.complement()
def transcribe(dna_seq):
        dna = Seq(dna_seq, generic_dna)

        #Fita de RNA:
        rna = dna.complement().transcribe()

def translate(seq_dna):
        #Fita de Dna:
        dna = Seq(seq_dna, generic_dna)

        #Fita de RNA mensageiro
        rna = dna.complement().transcribe()

        prot = rna.translate(to_stop=False, table=1)

#Copy functions
def copy_dnam():
        pyperclip.copy(entry_complementary_dna.get())

def copy_rna():
        pyperclip.copy(entry_rna.get())

def copy_seq_aa():
        pyperclip.copy(entry_seq_aa.get())
#functions Pairwise




#functions matrix pairwise



#
#buttons
button_run_dna = tk.Button(text='Run!', command=exit, bg='blue', fg='white', font=('helvetica', 9, 'bold'))
canvas1.create_window(400, 140, window=button_run_dna)#

button_copy_complementary_dna = tk.Button(text='Copy!', command=copy_dnam, bg='white', fg='black', font=('helvetica', 9, 'bold'))
canvas1.create_window(400, 180, window=button_copy_complementary_dna)#copy dna complementary

button_copy_rna = tk.Button(text='Copy!', command=copy_rna, bg='white', fg='black', font=('helvetica', 9, 'bold'))
canvas1.create_window(400, 220, window=button_copy_rna)#copy rna

button_copy_seq_aa = tk.Button(text='Copy!', command=copy_seq_aa, bg='white', fg='black', font=('helvetica', 9, 'bold'))
canvas1.create_window(400, 260, window=button_copy_seq_aa)# copy seq AA
#working until here
button_run_pairwisealignment = tk.Button(text='Run!', command=exit, bg = 'green', fg='white', font=('helvedica',9,'bold'))
canvas1.create_window(240,555, window=button_run_pairwisealignment)

button_run_matrix_alignment = tk.Button(text='Run!', command=exit, bg = 'green', fg='white', font=('helvedica',9,'bold'))
canvas1.create_window(530,460, window=button_run_matrix_alignment)



root.mainloop()
'''def copy():
    temp=inputBox.get()
    output.insert(END, temp)
    inputBox.delete(0, END)

#input box
inputBox = Entry(main, width=20, bg='white')
inputBox.grid(row=0, column=0, sticky=W)
#button
Button(main, width=20, text='copy', command=copy).grid(row=1, column=0, sticky=W)
#output box
output = Text(main, width=20, height=5, background='white')
output.grid(row=2,column=0,sticky=W)
main.mainloop()'''