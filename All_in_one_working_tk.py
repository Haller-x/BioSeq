####DNA COMPLEMENTAR###
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna# i need to change it
from Bio.SubsMat import MatrixInfo as matlist
from tkinter import filedialog
import tkinter.messagebox
import tkinter as tk
import pyperclip
from tkinter import ttk
from dict_values import *

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

####need to change it
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

entry_prot_seq1_matrix = tk.Entry (root)#seq aa
canvas1.create_window(485, 370, window=entry_prot_seq1_matrix)


entry_prot_seq2_matrix = tk.Entry (root)#seq aa
canvas1.create_window(485, 400, window=entry_prot_seq2_matrix)



def error_dna_codon():
        tkinter.messagebox.showerror('Warning', 'Sequence not a multiple of three')

def wrong_input():
        tkinter.messagebox.showerror('Error', 'Not a valid sequence')
        #maybe add some conference for ATGC
def generateTTR():
        try:
                dnaC = generate_dnc_c(entry_dna_seq.get())
                if int(len(entry_dna_seq.get())) % 3 != 0:
                        error_dna_codon()
                rna = generate_rna(entry_dna_seq.get())
                prot = generate_prot(entry_dna_seq.get())
                clear()
                entry_complementary_dna.insert(0,str(dnaC))
                entry_rna.insert(0,rna)
                entry_seq_aa.insert(0,prot)
        except:
                wrong_input()
def clear():
        entry_complementary_dna.delete(0, 'end')
        entry_rna.delete(0, 'end')
        entry_seq_aa.delete(0, 'end')

def generate_dnc_c(dna_seq):
        dnaM = Seq(dna_seq, generic_dna)
        dnaC = dnaM.complement()
        return dnaC

def generate_rna(dna_seq):
        dnaM = Seq(dna_seq, generic_dna)
        dnaC = dnaM.complement()
        rna = dnaM.complement().transcribe()
        return rna

def generate_prot(dna_seq):
        dnaM = Seq(dna_seq, generic_dna)
        dnaC = dnaM.complement()
        rna = dnaM.complement().transcribe()
        prot = rna.translate(to_stop=False, table=1)
        return prot

#Copy functions
def copy_dnam():
        pyperclip.copy(entry_complementary_dna.get())

def copy_rna():
        pyperclip.copy(entry_rna.get())

def copy_seq_aa():
        pyperclip.copy(entry_seq_aa.get())

#save funcions (bioseq)
def save_fasta_dnac():
        save_dnac = tk.Tk()
        save_dnac.withdraw()
        save_dnac.filename = filedialog.asksaveasfilename(title = "Save as ", defaultextension='.txt', filetypes = (("txt files","*.txt"),))
        with open(save_dnac.filename, "w") as save_dna_out:
                save_dna_out.write(entry_complementary_dna.get())

def save_fasta_rna():
        save_rna = tk.Tk()
        save_rna.withdraw()
        save_rna.filename = filedialog.asksaveasfilename(title = "Save as ", defaultextension='.txt', filetypes = (("txt files","*.txt"),))
        with open(save_rna.filename, "w") as save_rna_out:
                save_rna_out.write(entry_rna.get())

def save_fasta_seq_aa():
        save_seq_aa = tk.Tk()
        save_seq_aa.withdraw()
        save_seq_aa.filename = filedialog.asksaveasfilename(title = "Save as ", defaultextension='.txt', filetypes = (("txt files","*.txt"),))
        with open(save_seq_aa.filename, "w") as save_seq_aa_out:
                save_seq_aa_out.write(entry_seq_aa.get())


def pairwise_seq():
        try:
                seq1 = entry_seq_dna1_pairwise.get()
                seq2 = entry_seq_dna2_pairwise.get()
                match_pesos= entry_match_weight_pairwise.get()
                mismatch_pesos= entry_mismatch_weight_pairwise.get()
                gap_abertura= entry_gap_open.get()
                gap_extensao= entry_gap_extend.get()
                align= pairwise2.align.globalms(seq1,seq2,match_pesos,mismatch_pesos,gap_abertura,gap_extensao)
                save_pairwise = tk.Tk()
                save_pairwise.withdraw()
                save_pairwise.filename = filedialog.asksaveasfilename(title = "Save as ", defaultextension='.txt', filetypes = (("txt files","*.txt"),))
                with open(save_pairwise.filename, "w") as save_pairwise_out:
                        save_pairwise_out.write(pairwise2.format_alignment(*align[0]))
                return (pairwise2.format_alignment(*align[0]))
        except:
                tkinter.messagebox.showerror('Error', 'Not a valid sequence or values\nUse dot instead of comma')


def clear_pairwise():
        entry_seq_dna1_pairwise.delete(0, 'end')
        entry_seq_dna2_pairwise.delete(0, 'end')
        entry_match_weight_pairwise.delete(0, 'end')
        entry_mismatch_weight_pairwise.delete(0, 'end')
        entry_gap_open.delete(0, 'end')
        entry_gap_extend.delete(0, 'end')

def pairwise_default_seq(match_pesos=5,mismatch_pesos=-4,gap_abertura=-2,gap_extensao=-0.5):
        try:
                seq1 = entry_seq_dna1_pairwise.get()
                seq2 = entry_seq_dna2_pairwise.get()
                align= pairwise2.align.globalms(seq1,seq2,match_pesos,mismatch_pesos,gap_abertura,gap_extensao)

                save_pairwise_def = tk.Tk()
                save_pairwise_def.withdraw()
                save_pairwise_def.filename = filedialog.asksaveasfilename(title = "Save as ", defaultextension='.txt', filetypes = (("txt files","*.txt"),))
                with open(save_pairwise_def.filename, "w") as save_pairwise_def_out:
                        save_pairwise_def_out.write(pairwise2.format_alignment(*align[0]))
                return (pairwise2.format_alignment(*align[0]))
        except:
                tkinter.messagebox.showerror('Error', 'Not a valid sequence')

def clear_matrix_default():
        entry_prot_seq1_matrix.delete(0, 'end')
        entry_prot_seq2_matrix.delete(0, 'end')



#functions matrix pairwise


def alignment_matrix():
        try:
                seq1 = entry_prot_seq1_matrix.get()
                seq2 = entry_prot_seq2_matrix.get()
                chosen = dict_matrix[chosen_matrix]
                alinhamento =  pairwise2.align.globaldx(seq1,seq2,chosen)
                save_pairwise_prot = tk.Tk()
                save_pairwise_prot.withdraw()
                save_pairwise_prot.filename = filedialog.asksaveasfilename(title = "Save as ", defaultextension='.txt', filetypes = (("txt files","*.txt"),))
                with open(save_pairwise_prot.filename, "w") as save_pairwise_prot_out:
                        save_pairwise_prot_out.write(pairwise2.format_alignment(*alinhamento[0]))

        except:
                tkinter.messagebox.showerror('Error', 'Not a valid sequence')

#
#buttons
button_run_dna = tk.Button(text='Run!', command=generateTTR, bg='blue', fg='white', font=('helvetica', 9, 'bold'))
canvas1.create_window(400, 140, window=button_run_dna)#

button_copy_complementary_dna = tk.Button(text='Copy!', command=copy_dnam, bg='white', fg='black', font=('helvetica', 9, 'bold'))
canvas1.create_window(400, 180, window=button_copy_complementary_dna)#copy dna complementary

button_copy_rna = tk.Button(text='Copy!', command=copy_rna, bg='white', fg='black', font=('helvetica', 9, 'bold'))
canvas1.create_window(400, 220, window=button_copy_rna)#copy rna

button_copy_seq_aa = tk.Button(text='Copy!', command=copy_seq_aa, bg='white', fg='black', font=('helvetica', 9, 'bold'))
canvas1.create_window(400, 260, window=button_copy_seq_aa)# copy seq AA
##save buttons
button_save_dnaC = tk.Button(text='Save as txt!', command=save_fasta_dnac, bg='white', fg='black', font=('helvetica', 9, 'bold'))
canvas1.create_window(480, 180, window=button_save_dnaC)# copy seq AA

bbutton_save_rna = tk.Button(text='Save as txt!', command=save_fasta_rna, bg='white', fg='black', font=('helvetica', 9, 'bold'))
canvas1.create_window(480, 220, window=bbutton_save_rna)# copy seq AA

button_save_seq_aa = tk.Button(text='Save as txt!', command=save_fasta_seq_aa, bg='white', fg='black', font=('helvetica', 9, 'bold'))
canvas1.create_window(480,260, window=button_save_seq_aa)# copy seq AA

#run and clear buttons
button_run_pairwisealignment = tk.Button(text='Run!', command=pairwise_seq, bg = 'green', fg='white', font=('helvedica',9,'bold'))
canvas1.create_window(240,555, window=button_run_pairwisealignment)

button_run_pairwisealignment_default = tk.Button(text='Run default!', command=pairwise_default_seq, bg = 'green', fg='white', font=('helvedica',9,'bold'))
canvas1.create_window(180,555, window=button_run_pairwisealignment_default)

button_clear_pairwisealignment_default = tk.Button(text='Clear fields', command=clear_pairwise, bg = 'green', fg='white', font=('helvedica',9,'bold'))
canvas1.create_window(100,555, window=button_clear_pairwisealignment_default)

button_run_matrix_alignment = tk.Button(text='Run!', command=alignment_matrix, bg = 'green', fg='white', font=('helvedica',9,'bold'))
canvas1.create_window(530,460, window=button_run_matrix_alignment)

button_clear_matrix_alignment = tk.Button(text='Clear fields', command=clear_matrix_default, bg = 'green', fg='white', font=('helvedica',9,'bold'))
canvas1.create_window(470,460, window=button_clear_matrix_alignment)



root.mainloop()
