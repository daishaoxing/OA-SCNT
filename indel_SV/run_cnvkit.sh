ref="/data/dsx/Genome/fasta_gtf/Macaca_mulatta.Mmul_8.0.1.dna.chr.fa"
refFlat="/data/dsx/Genome/fasta_gtf/Macaca_mulatta.Mmul_8.0.1.96.onlychr.refFlat"
cnvkit.py batch *.bam -n -m wgs -f $ref --annotate $refFlat --scatter --diagram -d cnvkit_outputv2 -p 40
cnvkit.py heatmap Donor-cell1.cns \
				Donor-cell2.cns \
				Donor-cell3.cns \
				Donor-cell4.cns \
				SCNT1.cns \
				SCNT2.cns \
				SCNT3.cns \
				SCNT4.cns \
				SCNT5.cns \
				SCNT6.cns \
				SCNT7.cns \
				SCNT8.cns \
				SCNT9.cns \
				SCNT10-Body.cns \
				SCNT11-Head.cns \
				SCNT-ABE1.cns \
				SCNT-ABE2.cns \
				SCNT-ABE3.cns \
				SCNT-ABE4.cns \
				SCNT-ABE5.cns \
				SCNT-ABE6.cns \
				SCNT-ABE7.cns \
				SCNT-ABE8.cns \
				SCNT-ABE9-Muscle.cns \
				SCNT-ABE10-Skin.cns \
				SCNT-ABE11-Stomach.cns \
				-o CNVheatmap.pdf -x Female
