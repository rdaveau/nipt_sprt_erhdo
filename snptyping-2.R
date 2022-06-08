#!/usr/bin/Rscript
rm(list = ls())
options(warn = -1, bitmapType = 'cairo')
X11.options(type = 'cairo')

rlibs = c('optparse', 'RColorBrewer', 'ggplot2', 'reshape2', 'stringr')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

proband = list('AC' = 'Affected child', 'UC' = 'Unaffected child', 'AR' = 'Affected close relative', 'UR' = 'Unaffected close relative')

opts = parse_args(object = (args = OptionParser(option_list = list(
	make_option(opt_str = '--pref', type = 'character', help = 'Input files prefix name'),
	make_option(opt_str = '--idir', type = 'character', help = 'Path to input directory'),
	make_option(opt_str = '--odir', type = 'character', help = 'Path to output directory'),
	make_option(opt_str = '--bed', type = 'character', help = 'Path to design BED file'),
	make_option(opt_str = '--genmode', type = 'character', help = 'Genetic inheritance (ADMI, ADPI or ARI)'),
	make_option(opt_str = '--proband', type = 'character', help = 'Proband type (AC, UC, AR or UR)'),
	make_option(opt_str = '--matmut', type = 'character', help = 'Maternal mutation (if any)'),
	make_option(opt_str = '--patmut', type = 'character', help = 'Paternal mutation (if any)'),
	make_option(opt_str = '--HTZ', type = 'double', default = .15, help = 'Minimal allele ratio to call an heterozygote site'),
	make_option(opt_str = '--min.NDP', type = 'integer', default = 8, help = 'Minimal sequencing depth in nuclear dna'),
	make_option(opt_str = '--min.PDP', type = 'integer', default = 15, help = 'Minimal sequencing depth in plasmatic dna'),
	make_option(opt_str = '--min.SNP', type = 'integer', default = 10, help = 'Minimal number of SNP3 or -4 for haplotype inference estimation'),
	make_option(opt_str = '--min.dist', type = 'integer', default = 50, help = 'Minimal distance between 2 contiguous SNPs for SPRT analysis'),
	make_option(opt_str = '--alpha', type = 'double', default = .001, help = 'Alpha risk for SPRT analysis'),
	make_option(opt_str = '--beta', type = 'double', default = .001, help = 'Beta risk for SPRT analysis'),
	make_option(opt_str = '--ZFXY', type = 'character', help = 'Path to ZFX/ZFY coverage file in maternal plasmatic DNA')
), usage = 'snptyping.R [options]')))

min.depth = paste0(paste(opts$min.NDP, opts$min.PDP, sep = '/'), 'X')
snpmap = read.table(file = file.path(opts$idir, paste0(opts$pref, '_snpmap.tsv')), sep = "\t", header = TRUE, as.is = FALSE)
bedfile = read.table(file = opts$bed, sep = "\t", header = FALSE, as.is = TRUE, col.names = c('CHR', 'START', 'END', 'NAME'))
ifiles = list.files(path = opts$idir, pattern = paste0(opts$pref, '_chr'))
chr = gsub(pattern = paste0(opts$pref, "_chr(.+)\\.tsv$"), replacement = "\\1", x = ifiles, perl = TRUE)

SNP.TABLE = matrix(data = 0, nrow = length(levels(snpmap$SNP.TYPE)), ncol = length(levels(snpmap$SNP.SUBTYPE)))
SNP1 = c()
MAF = c()
ERROR.RATE = list()
SNP3 = list()
SNP4 = list()

# Extract mat-/paternal mutation coordinates and target gene(s) (if any)
mut = sapply(c('matmut', 'patmut'), function(mut){
	o = list()
	if(opts[[mut]] != 'NA'){
		x = unlist(strsplit(x = opts[[mut]], split = ':'))
		o$chr = x[1]
		o$range = as.numeric(unlist(str_extract_all(string = x[2], pattern = "\\d+")))
		if(length(o$range) != 2){
			o$range = rep(o$range, 2)
		}
		m = subset(bedfile, CHR == o$chr)
		o$target = unique(m[c(rev(which(o$range[1] >= m$START))[1] : which(o$range[2] <= m$END)[1]), 'NAME'])
	}
	o
}, simplify = FALSE)

DPNI_DP = c() # Mean coverage in maternal plasmatic DNA

for(ii in 1:length(chr)){

	dat = read.table(file = file.path(opts$idir, ifiles[ii]), sep = "\t", header = TRUE, as.is = TRUE)

	if(nrow(dat)){
		DPNI_DP = c(DPNI_DP, mean(dat$DPNI.TRC)) # Avoid NaN in mean coverage calculation when zero-SNP
	}

	for(x in c('SNP.TYPE', 'SNP.SUBTYPE')){
		dat[[x]] = factor(dat[[x]], levels = levels(snpmap[[x]]))
	}

	# SNP types distribution
	SNP.TABLE = SNP.TABLE + (m = with(dat, table(SNP.TYPE, SNP.SUBTYPE)))
	n = nrow(m)
	m = 100 * m / rowSums(m)
	m[!is.finite(m)] = 0

	m = melt(m)
	p = ggplot(data = m) +
		geom_bar(aes(x = SNP.TYPE, y = value, fill = SNP.SUBTYPE), stat = 'identity') +
		annotate(geom = 'text', x = 1:n, y = 1, label = paste0(table(dat$SNP.TYPE), "\n")) +
		scale_fill_manual(name = 'subtype', values = c(brewer.pal(9, 'Set1'), brewer.pal(8, 'Set2')[7:1])) +
		labs(
			title = paste(opts$pref, 'distribution of SNP', min.depth),
			subtitle = paste(paste0('chr', chr[ii]), opts$genmode, proband[[opts$proband]], paste0('(n=', nrow(dat), ')')),
			x = 'SNP types',
			y = 'SNP counts (%)'
		) + theme_minimal() + theme(
			plot.title = element_text(size = rel(1.8), lineheight = 2, vjust = 1, face = 'bold'), plot.subtitle = element_text(size = rel(1.5)),
			legend.title = element_text(size = rel(1.5)), legend.position = 'right', legend.direction = 'vertical',
			axis.title = element_text(size = rel(1.5)), axis.text = element_text(size = rel(1.5))
		)
	ggsave(filename = file.path(opts$odir, paste(opts$pref, paste0('chr', chr[ii]), 'SNP.TABLE.png', sep = '_')), plot = p, width = 10, height = 8, device = 'png', dpi = 150)

	# SNP1-based ff estimates
	m = subset(dat, SNP.TYPE == 1)
	if(nrow(m) >= opts$min.SNP){
		x = m$DPNI.RRC
		i = m$Mo.GT == 'AA'
		x[i] = m$DPNI.ARC[i]
		x = 100 * x / m$DPNI.TRC
		SNP1 = c(SNP1, x)
		dx = density(x)
		xx = dx$x[which.max(dx$y)]
		p = ggplot(data = melt(x)) +
			geom_density(aes(value), fill = brewer.pal(9, 'Set1')[2], colour = brewer.pal(9, 'Set1')[2], alpha = .1) +
			geom_vline(xintercept = c(dx$x[which.max(dx$y)] - sd(x), dx$x[which.max(dx$y)], dx$x[which.max(dx$y)] + sd(x)), linetype = c('dotted', 'dashed', 'dotted'), alpha = .75) +
			geom_label(data = melt(xx), aes(x = value, y = 0), label = paste0('ff=', round(2 * xx, digits = 1), '%'), size = 5) +
			labs(
				title = paste(opts$pref, 'ff density estimates', min.depth),
				subtitle = paste(paste0('chr', chr[ii]), opts$genmode, proband[[opts$proband]], paste0('(#SNP1=', length(x), ')')),
				x = 'Fetal fraction (%)',
				y = 'Density estimates'
			) + theme_minimal() + theme(
				plot.title = element_text(size = rel(1.8), lineheight = 2, vjust = 1, face = 'bold'), plot.subtitle = element_text(size = rel(1.5)),
				legend.title = element_text(size = rel(1.5)), legend.position = 'right', legend.direction = 'vertical',
				axis.title = element_text(size = rel(1.5)), axis.text = element_text(size = rel(1.5))
			)
		ggsave(filename = file.path(opts$odir, paste(opts$pref, paste0('chr', chr[ii]), 'FF_PAF_ESTIMATES.png', sep = '_')), plot = p, width = 10, height = 8, device = 'png', dpi = 150)
	}

	# MAF-based ff estimates
	i = grep('W', dat$SNP.TYPE, invert = TRUE)
	MAF = c(MAF, 100 * apply(subset(dat[i, ], select = c(DPNI.RRC, DPNI.ARC)), 1, min) / dat[i, ]$DPNI.TRC)
	MAF = MAF[!!MAF]

	# SNP2-based sequencing error rate estimates
	m = subset(dat, SNP.TYPE == 2)
	if(nrow(m)){
		ERROR.RATE[[paste0('chr', chr[[ii]])]] = data.frame(sapply(c('Mo', 'Fa', 'CI', 'DPN', 'DPNI'), function(x){
			100 * m[[paste0(x, '.BRC')]] / m[[paste0(x, '.TRC')]]
		}))
	}

	# SNP3-based haplotype inference
	m = subset(dat, SNP.TYPE == 3)
	mm = subset(bedfile, CHR == chr[ii])
	m$BED.TARGET = sapply(m$POS, function(x){
		unique(mm$NAME[x >= mm$START & x <= mm$END])
	})
	if(nrow(m) >= opts$min.SNP){
		SNP3[[paste0('chr', chr[[ii]])]] = m
	}

	# SNP4-based haplotype inference
	m = dat[grep(pattern = '4', x = dat$SNP.TYPE), ]
	mm = subset(bedfile, CHR == chr[ii])
	m$BED.TARGET = sapply(m$POS, function(x){
		unique(mm$NAME[x >= mm$START & x <= mm$END])
	})
	if(nrow(m) >= opts$min.SNP){
		m$SNP.TYPE.BAK = m$SNP.TYPE
		SNP4[[paste0('chr', chr[[ii]])]] = m
	}
}

# SNP types distribution
n = nrow(SNP.TABLE)
m = 100 * SNP.TABLE / rowSums(SNP.TABLE)
m[!is.finite(m)] = 0
m = melt(m)
p = ggplot(data = m) +
	geom_bar(aes(x = SNP.TYPE, y = value, fill = SNP.SUBTYPE), stat = 'identity') +
	annotate(geom = 'text', x = 1:n, y = 1, label = paste0(rowSums(SNP.TABLE), "\n")) +
	scale_fill_manual(name = 'subtype', values = c(brewer.pal(9, 'Set1'), brewer.pal(8, 'Set2')[7:1])) +
	labs(
		title = paste(opts$pref, 'distribution of SNP', min.depth),
		subtitle = paste(opts$genmode, proband[[opts$proband]], paste0('(n=', sum(SNP.TABLE), ')')),
		x = 'SNP types',
		y = 'SNP counts (%)'
	) + theme_minimal() + theme(
		plot.title = element_text(size = rel(1.8), lineheight = 2, vjust = 1, face = 'bold'), plot.subtitle = element_text(size = rel(1.5)),
		legend.title = element_text(size = rel(1.5)), legend.position = 'right', legend.direction = 'vertical',
		axis.title = element_text(size = rel(1.5)), axis.text = element_text(size = rel(1.5))
	)
ggsave(filename = file.path(opts$odir, paste(opts$pref, 'SNP.TABLE.png', sep = '_')), plot = p, width = 10, height = 8, device = 'png', dpi = 150)
write.table(x = cbind(SNP.TYPE = rownames(SNP.TABLE), SNP.TABLE), file = file.path(opts$odir, paste(opts$pref, 'SNP.TABLE.tsv', sep = '_')), sep = "\t", quote = FALSE, row.names = FALSE)

# SNP1-based ff estimates
FF.TABLE = data.frame('PAF' = NA)
if(!is.null(SNP1)){
	dx = density(SNP1)
	ff = xx = dx$x[which.max(dx$y[dx$x < 25])]

	FF.TABLE$PAF = format(2 * xx, digits = 1, nsmall = 1, trim = TRUE)

	p = ggplot() +
		geom_density(data = melt(SNP1), aes(value), fill = brewer.pal(9, 'Set1')[2], colour = brewer.pal(9, 'Set1')[2], alpha = .1) +
		geom_vline(xintercept = c(xx - sd(SNP1), xx, xx + sd(SNP1)), linetype = c('dotted', 'dashed', 'dotted'), alpha = .75) +
		geom_label(data = melt(xx), aes(x = value, y = 0), label = paste0('ff=', round(2 * xx, digits = 1), '%'), size = 5) +
		labs(
			title = paste(opts$pref, 'ff density estimates', min.depth),
			subtitle = paste(opts$genmode, proband[[opts$proband]], paste0('(#SNP1=', length(SNP1), ')')),
			x = 'Paternal Allele Frequency (%)',
			y = 'Density estimates'
		) + theme_minimal() + theme(
			plot.title = element_text(size = rel(1.8), lineheight = 2, vjust = 1, face = 'bold'), plot.subtitle = element_text(size = rel(1.5)),
			legend.title = element_text(size = rel(1.5)), legend.position = 'right', legend.direction = 'vertical',
			axis.title = element_text(size = rel(1.5)), axis.text = element_text(size = rel(1.5))
		)
	ggsave(filename = file.path(opts$odir, paste(opts$pref, 'FF_PAF_ESTIMATES.png', sep = '_')), plot = p, width = 10, height = 8, device = 'png', dpi = 150)
}

# MAF-based ff estimates
dx = density(MAF)
xx = dx$x[which.max(dx$y[dx$x < 25])]

FF.TABLE$MAF = format(2 * xx, digits = 1, nsmall = 1, trim = TRUE)

# Add FF estimate from ZFX/ZFY coverage
covfile = read.table(file = opts$ZFXY, sep = "\t", header = FALSE, as.is = TRUE, col.names = c('CHR', 'START', 'END', 'NAME', 'POS', 'DP'))
ZFXY = list(mu = sapply(split(x = covfile, f = covfile$NAME), function(x){ mean(x$DP) }))
ZFXY$ff = ifelse(test = ZFXY$mu['ZFY'] < 1, yes = NA, no = 200 * ZFXY$mu['ZFY'] / sum(ZFXY$mu))
FF.TABLE$ZFXY = ifelse(test = is.na(ZFXY$ff), yes = NA, no = format(ZFXY$ff, digits = 1, nsmall = 1, trim = TRUE))

write.table(x = FF.TABLE, file = file.path(opts$odir, paste(opts$pref, 'FF.TABLE.tsv', sep = '_')), sep = "\t", quote = FALSE, row.names = FALSE)

# Write ZFY/ZFY and PLASMA coverage data
write.table(x = lapply(c(ZFXY$mu, PLASMA = mean(DPNI_DP)), function(x) format(x = x, digits = 1, nsmall = 1, trim = TRUE)),
	file = file.path(opts$odir, paste(opts$pref, 'COV.TABLE.tsv', sep = '_')), sep = "\t", quote = FALSE, row.names = FALSE)

# Set ff based on MAF only
ff = 2 * xx
ff.sd = sd(MAF[MAF < 100*opts$HTZ])

p = ggplot() +
	geom_density(data = melt(MAF), aes(value), fill = brewer.pal(9, 'Set1')[2], colour = brewer.pal(9, 'Set1')[2], alpha = .1) +
	geom_vline(xintercept = c(xx - sd(MAF[MAF < 100*opts$HTZ]), xx, xx + sd(MAF[MAF < 100*opts$HTZ])), linetype = c('dotted', 'dashed', 'dotted'), alpha = .75) +
	geom_label(data = melt(xx), aes(x = value, y = 0), label = paste0('ff=', round(2 * xx, digits = 1), '%'), size = 5) +
	labs(
		title = paste(opts$pref, 'ff density estimates', min.depth),
		subtitle = paste(opts$genmode, proband[[opts$proband]], paste0('(#SNP=', length(MAF), ')')),
		x = 'Minor Allele Frequency (%)',
		y = 'Density estimates'
	) + theme_minimal() + theme(
		plot.title = element_text(size = rel(1.8), lineheight = 2, vjust = 1, face = 'bold'), plot.subtitle = element_text(size = rel(1.5)),
		legend.title = element_text(size = rel(1.5)), legend.position = 'right', legend.direction = 'vertical',
		axis.title = element_text(size = rel(1.5)), axis.text = element_text(size = rel(1.5))
	)
ggsave(filename = file.path(opts$odir, paste(opts$pref, 'FF_MAF_ESTIMATES.png', sep = '_')), plot = p, width = 10, height = 8, device = 'png', dpi = 150)

# SNP2-based sequencing error rate estimates
ERROR.RATE = do.call(rbind, ERROR.RATE)
names(ERROR.RATE)[grep(pattern = 'DPNI', names(ERROR.RATE))] = 'NIPD'
m = melt(ERROR.RATE)
dx = density(m$value[(i = m$value!=0 & m$value<5)])
xx = dx$x[which.max(dx$y)]
p = ggplot() +
	geom_density(data = m[i, ], aes(value, fill = variable, colour = variable), alpha = .1) +
	geom_density(data = m[i, ], aes(value), fill = NA, colour = brewer.pal(9, 'Set1')[1], alpha = .1, linetype = 'dashed') +
	geom_vline(xintercept = c(xx - sd(m$value), xx, xx + sd(m$value)), linetype = c('dotted', 'dashed', 'dotted'), alpha = .75) +
	geom_label(data = melt(xx), aes(x = value, y = 0), label = paste0('Error=', round(xx, digits = 1), '%'), size = 5) +
	labs(
		title = paste(opts$pref, 'Error rate density estimates', min.depth),
		subtitle = paste(opts$genmode, proband[[opts$proband]], paste0('(#SNP2=', nrow(ERROR.RATE), ')')),
		x = 'Error rate (%)',
		y = 'Density estimates',
		fill = '', colour = ''
	) + theme_minimal() + theme(
		plot.title = element_text(size = rel(1.8), lineheight = 2, vjust = 1, face = 'bold'), plot.subtitle = element_text(size = rel(1.5)),
		legend.title = element_text(size = rel(1.5)), legend.position = 'right', legend.direction = 'vertical',
		axis.title = element_text(size = rel(1.5)), axis.text = element_text(size = rel(1.5))
	) + xlim(0, 5)
ggsave(filename = file.path(opts$odir, paste(opts$pref, 'ERROR_RATE_ESTIMATES.png', sep = '_')), plot = p, width = 10, height = 8, device = 'png', dpi = 150)

# SNP3-based haplotype inference
min.VAF = max(3*xx, .5*ff-ff.sd)

if(opts$genmode != 'ADMI'){
	invisible(lapply(names(SNP3), function(chr){

		x = SNP3[[chr]]
		j = setdiff(1:nrow(x), (i = which(x$Fa.ALLELE == 'A')))

		x$DPNI.VAF = x$DPNI.HAP = NA
		x$DPNI.VAF[i] = 100 * x$DPNI.RRC[i] / x$DPNI.TRC[i]
		x$DPNI.VAF[j] = 100 * x$DPNI.ARC[j] / x$DPNI.TRC[j]
		jj = setdiff(1:nrow(x), (ii = which(x$DPNI.VAF >= min.VAF)))
		x$DPNI.HAP[ii] = x$Fa.HAP[ii]
		x$DPNI.HAP[jj] = 0

		x$DPN.HAP = sapply(1:nrow(x), function(k){
			ifelse(length(grep(pattern = x$Fa.ALLELE[k], x = x$DPN.GT[k])), x$Fa.HAP[k], ifelse(x$Fa.HAP[k] == 3, 4, 3))
		})

		x = split(x = x, f = x$BED.TARGET)

		lapply(setdiff(names(x), 'Ethnie'), function(NAME){

			max.DPNI.VAF = max(x[[NAME]]$DPNI.VAF)
			xx = split(x = x[[NAME]], f = x[[NAME]]$SNP.SUBTYPE)

			# Merge SNP3 relatives ie. A+D+E+H and B+C+F+G
			xx = lapply(list(
				'A' = do.call(rbind, list(xx$A, xx$D, xx$E, xx$H)),
				'B' = do.call(rbind, list(xx$B, xx$C, xx$F, xx$G))
			), function(x) x[order(x$POS), ])

			xx.label = list(A = 'A/D/E/H', B = 'B/C/F/G')

			lapply(names(xx), function(SUBTYPE){

				m = subset(xx[[SUBTYPE]], select = c(POS, DPN.HAP, DPNI.HAP, DPNI.VAF))
				if(nrow(m) >= opts$min.SNP){
					m$DPNI.HAP = factor(m$DPNI.HAP, levels = c(0, 3:4), labels = c(paste('?', paste0('(', sum(!m$DPNI.HAP), ')')), paste('III', paste0('(', sum(m$DPNI.HAP == 3), ')')), paste('IV', paste0('(', sum(m$DPNI.HAP == 4), ')'))))
					m$DPN.HAP = factor(m$DPN.HAP, levels = c(0, 3:4), labels = levels(m$DPNI.HAP))
					p = ggplot(data = m) +
						geom_tile(aes(x = 1:nrow(m), y = DPNI.VAF / max.DPNI.VAF, fill = DPNI.HAP), height = .1) +
						geom_tile(aes(x = 1:nrow(m), y = -.2, fill = DPN.HAP), height = .15) +
						geom_line(data = m, aes(x = 1:nrow(m), y = scales::rescale(x = c(0, diff(POS)), to = c(-0.12, 1)))) +
						geom_hline(yintercept = min.VAF / max.DPNI.VAF, linetype = 'dashed', alpha = .75) +
						geom_text(data = data.frame(x = nrow(m), y = -.2, label = 'DPN'), aes(x = x, y = y, label = label), hjust = 0, vjust = .5, nudge_x = 1) +
						geom_label(data = melt(min.VAF / max.DPNI.VAF), aes(x = nrow(m), y = value), label = paste0('min=', format(min.VAF, digits = 2, nsmall = 2, trim = TRUE), '%'), size = 3) +
						scale_fill_manual(name = 'Inherited\nHAP', values = brewer.pal(9, 'Set1')[c(9,1,3)], drop = FALSE) +
						labs(
							title = paste(opts$pref, paste(chr, NAME, sep = ':'), 'Fa haplotype inheritance', min.depth),
							subtitle = paste(opts$genmode, proband[[opts$proband]], paste0('(#SNP3', xx.label[[SUBTYPE]], '=', nrow(m), ')')),
							x = '',
							y = 'VAF (AU)'
						) + theme_minimal() + theme(
							plot.title = element_text(size = rel(1.8), lineheight = 2, vjust = 1, face = 'bold'), plot.subtitle = element_text(size = rel(1.1)),
							legend.title = element_text(size = rel(1.5)), legend.position = 'right', legend.direction = 'vertical',
							axis.title = element_text(size = rel(1.5)), axis.text = element_blank()
						) + ylim(-.4, 1)

					# Add genomic coordinates
					xt = c(1, c(order(rev(rev(diff(m$POS))[-1]), decreasing = TRUE) + 1)[1:5], nrow(m))
					p = p + geom_text(data = data.frame(x = xt, y = -.125, z = m$POS[xt]), aes(x = x, y = y, label = z), size = 2.5, angle = 45, hjust = 1, vjust = .5)

					if(!!length(mut$patmut)){
						if(paste0('chr', mut$patmut$chr) == chr & !is.na(match(NAME, mut$patmut$target))){
							i = c(rev(which(mut$patmut$range[1] >= m$POS))[1], which(mut$patmut$range[2] <= m$POS)[1])
							if(is.na(i[1])){
								i[1] = 1
							}
							if(is.na(i[2])){
								i[2] = nrow(m)
							}
							i = i - .5
							p = p +
								geom_segment(aes(x = i[1], xend = i[2], y = -.4, yend = -.4), size = 2) +
								geom_vline(xintercept = i, linetype = 'dashed')
						}
					}
					ggsave(filename = file.path(opts$odir, paste(opts$pref, chr, NAME, paste0('SNP3', SUBTYPE), 'INHERITANCE.png', sep = '_')), plot = p, width = 12, height = 3, device = 'png', dpi = 150)

					xx[[SUBTYPE]][['DPNI.VAF']] = format(x = xx[[SUBTYPE]][['DPNI.VAF']], digits = 3, nsmall = 3, trim = TRUE)
					for(a in c('Fa', 'DPN', 'DPNI')){
						xx[[SUBTYPE]][[paste0(a, '.HAP')]] = factor(xx[[SUBTYPE]][[paste0(a, '.HAP')]], levels = c(0, 3:4), labels = c('?', 'III', 'IV'))
					}
					write.table(x = xx[[SUBTYPE]], file = file.path(opts$odir, paste(opts$pref, chr, NAME, paste0('SNP3', SUBTYPE), 'INHERITANCE.tsv', sep = '_')), sep = "\t", quote = FALSE, row.names = FALSE)
				}
			})
		})
	}))
}

# SNP4-based haplotype inference (sprt)
if(opts$genmode != 'ADPI'){

	f = ff / 100

	q = list(
		q0 = list('4a' = .5 , '4b' = .5 * (1-f)),
		q1 = list('4a' = .5 * (1+f), '4b' = .5)
	)

	q_bounds = function(dp, q0, q1){
		d = (1-q1) / (1-q0)
		g = (q1 * (1-q0)) / (q0 * (1-q1))
		lapply(list('LOWER.QBOUND' = opts$beta / (1 - opts$alpha), 'UPPER.QBOUND' = (1 - opts$beta) / opts$alpha), function(x){
			(log(x) / dp - log(d)) / log(g)
		})
	}

	SPRT.TABLE = do.call(rbind, lapply(names(SNP4), function(chr){

		x = SNP4[[chr]]
		x = split(x = x, f = x$BED.TARGET)

		SPRT.TABLE = do.call(rbind, lapply(setdiff(names(x), 'Ethnie'), function(NAME){

			m = x[[NAME]]
			m = m[c(TRUE, diff(m$POS) >= opts$min.dist), ]
			M = split(x = m, f = as.character(m$SNP.TYPE))

			m = do.call(rbind, lapply(M, function(m){

				m$RELIABLE.SNP = with(m, DPN.GT == Mo.GT | DPN.GT == Fa.GT)
				m$Fa.DPN.ALLELE = ifelse(m$Fa.GT == 'AA', 'A', 'B')
				m$Mo.DPN.ALLELE = ifelse(m$DPN.GT == 'AA', 'A', ifelse(m$DPN.GT == 'BB', 'B', ifelse(m$Fa.DPN.ALLELE == 'A', 'B', 'A')))
				m$DPN.HAP = ifelse(m$Mo.DPN.ALLELE == m$RISK.HAP, 1, 2)
				m$DPN.HAP[!m$RELIABLE.SNP] = 0
				m$DPNI.RISK.ALLELE = paste0(ifelse(m$RISK.HAP == 'B', 'A', 'R'), 'RC')
				p = list(FWD = 1:nrow(m), REV = nrow(m):1)
				for(dir in names(p)){
					pp = p[[dir]]
					for(x in c('.PROP.RISK.HAP', '.CUM.DEPTH', '.LOWER.QBOUND', '.UPPER.QBOUND')){
						m[[paste0(dir, x)]] = NA
					}
					for(x in c('HAP', 'HAP.SMOOTHED')){
						m[[paste('DPNI', dir, x, sep = '.')]] = 0
					}
					m[[paste('SPRT', dir, 'BLOCK', sep = '.')]] = NA
					block_id = 1
					n = pc = 0
					j = c()
					for(i in 1:length(pp)){
						j = c(j, pp[i])
						pc = m[[paste0(dir, '.PROP.RISK.HAP')]][pp[i]] = (pc * n + m[[paste0('DPNI.', m$DPNI.RISK.ALLELE[pp[i]])]][pp[i]]) / (n + m$DPNI.TRC[pp[i]])
						n = m[[paste0(dir, '.CUM.DEPTH')]][pp[i]] = n + m$DPNI.TRC[pp[i]]
						qb = q_bounds(n, q$q0[[as.character(m$SNP.TYPE[pp[i]])]], q$q1[[as.character(m$SNP.TYPE[pp[i]])]])
						for(x in names(qb)){
							m[[paste(dir, x, sep = '.')]][pp[i]] = qb[[x]]
						}
						m[[paste('SPRT', dir, 'BLOCK', sep = '.')]][pp[i]] = block_id
						m[[paste('DPNI', dir, 'HAP', sep = '.')]][pp[i]] = ifelse(pc <= qb$LOWER.QBOUND, 2, ifelse(pc >= qb$UPPER.QBOUND, 1, 0))
						if(m[[paste('DPNI', dir, 'HAP', sep = '.')]][pp[i]] != 0){
							m[[paste('DPNI', dir, 'HAP.SMOOTHED', sep = '.')]][j] = m[[paste('DPNI', dir, 'HAP', sep = '.')]][pp[i]]
							block_id = block_id + 1
							n = pc = 0
							j = c()
						}
					}
					m[[paste('SPRT', dir, 'BLOCK', sep = '.')]] = paste(m$SNP.TYPE, m[[paste('SPRT', dir, 'BLOCK', sep = '.')]], sep = '.')
				}
				m
			}))
			m = m[order(m$POS), ]

			# plot sprt profile

			# Compute %overall identity between DPNI and DPN HAP
			# Original 3-factors PC calculation (same as pc in snptyping-1.R)
			pc_3f = sapply(c('FWD', 'REV'), function(x){

				i = apply((xx = m[, c('DPN.HAP', paste('DPNI', x, 'HAP.SMOOTHED', sep = '.'))]), 1, min) != 0
				if(!sum(i)){
					list(pc = NA, nc = NA)
				}else{
					list(
						pc = 100 * sum(sapply(1:nrow((xx = apply(xx[i, ], 1:2, function(x) substr(x = x, start = 1, stop = 1)))), function(i) length(unique(unlist(xx[i, ]))) == 1)) / nrow(xx),
						nc = 100 * sum(!i) / nrow(m)
					)
				}
			}, simplify = FALSE)
			pc_3f$ALL = sapply(c('pc', 'nc'), function(x){
				mean(c(pc_3f$FWD[[x]], pc_3f$REV[[x]]))
			}, simplify = FALSE)

			# 2-factors PC ie. w/o taking into account the DPN sample
			i = m[, 'DPNI.REV.HAP.SMOOTHED'] != 0
			j = m[, 'DPNI.FWD.HAP.SMOOTHED'] != 0
			nc_2f = mean(c(100 * sum(!i) / nrow(m), 100 * sum(!j) / nrow(m)))
			i = apply((xx = m[, c('DPNI.REV.HAP.SMOOTHED', 'DPNI.FWD.HAP.SMOOTHED')]), 1, min) != 0
			# nc_2f = 100 * sum(!i) / nrow(m)
 			if(length(which(i == TRUE)) > 0){
      			pc_2f = 100 * sum(sapply(1:nrow((xx = apply(xx[i, ], 1:2, function(x) substr(x = x, start = 1, stop = 1)))), function(i) length(unique(unlist(xx[i, ]))) == 1)) / nrow(xx)
  			}else{
    			pc_2f = NA
  			}

			mm = subset(m, select = c(DPNI.FWD.HAP.SMOOTHED, DPNI.REV.HAP.SMOOTHED, DPN.HAP))

			names(mm) = c('FWD', 'REV', 'DPN')
			mm = data.frame(
				x = rep(1:nrow(mm), 3),
				y = factor(rep(names(mm), each = nrow(mm)), levels = c('DPN','REV','FWD')),
				z = factor(unlist(mm), levels = c(0:2), labels = c('?', 'I', 'II'))
			)

			# Number of SPRT blocks
			nb = melt((n = sapply(c('FWD', 'REV'), function(x){
				length(unique(m[[paste('SPRT', x, 'BLOCK', sep = '.')]]))
			}, simplify = FALSE)))

			# Mean #SNP per block and type
			ns = sapply(c('4a', '4b'), function(x){
				xx = subset(m, SNP.TYPE == x)
				ifelse(!nrow(xx), 0, mean(sapply(list(xx$SPRT.FWD.BLOCK, xx$SPRT.REV.BLOCK), function(x) mean(table(x)))))
			}, simplify = FALSE)

			# Add concordance and block score
			# Original 3-factors CS calculation (same as c_score in snptyping-1.R)
			cs_3f = pc_3f$ALL$pc
			cs_3f = 1 - logb(100 - cs_3f, base = 100)
			if(length(which(cs_3f > 1)) > 0){
				cs_3f[cs_3f > 1] = 1
			}
			cs_3f[is.na(cs_3f)] = 0
			cs_3f = cs_3f - pc_3f$ALL$nc / 100
			cs_3f[which(cs_3f < 0)] = 0

			# 2-factors CS ie. w/o taking into account the DPN sample
			cs_2f = pc_2f
			cs_2f = 1 - logb(100 - cs_2f, base = 100)
			if(length(which(cs_2f > 1)) > 0){
				cs_2f[cs_2f > 1] = 1
			}
			cs_2f[is.na(cs_2f)] = 0
			cs_2f = cs_2f - nc_2f / 100
			cs_2f[which(cs_2f < 0)] = 0

			b_score = min(1 - logb(x =  100 * mean(unlist(ns)) / nrow(m), base = 100), 1)

			# Add SPRT summary to SPRT.TABLE
			SPRT.TABLE = data.frame(
				'GENE' = NAME, # Gene Name
				'#SNP4' = nrow(m), # Total #SNP4
				'#FWD' = n$FWD, # Total #FWD blocks
				'#REV' = n$REV, # Total #REV blocks
				'#MUA' = format(x = ns[['4a']], digits = 1, nsmall = 1, trim = TRUE), # Mean #SNP4a / block
				'#MUB' = format(x = ns[['4b']], digits = 1, nsmall = 1, trim = TRUE), # Mean #SNP4b / block
				'%FWD' = format(x = pc_3f$FWD$pc, digits = 1, nsmall = 1, trim = TRUE),
				'%REV' = format(x = pc_3f$REV$pc, digits = 1, nsmall = 1, trim = TRUE),
				'PC2' = format(x = pc_2f, digits = 1, nsmall = 1, trim = TRUE),
				'NC2' = format(x = nc_2f, digits = 1, nsmall = 1, trim = TRUE),
				'PC3' = format(x = pc_3f$ALL$pc, digits = 1, nsmall = 1, trim = TRUE),
				'NC3' = format(x = pc_3f$ALL$nc, digits = 1, nsmall = 1, trim = TRUE),
				'CS2' = format(x = cs_2f, digits = 2, nsmall = 2, trim = TRUE),
				'CS3' = format(x = cs_3f, digits = 2, nsmall = 2, trim = TRUE),
				'BS' = format(x = b_score, digits = 2, nsmall = 2, trim = TRUE),
			check.names = FALSE)

			p = ggplot(data = mm, aes(x = x, y = y)) +
				geom_tile(aes(fill = z), height = .7) +

				scale_fill_manual(name = 'Inherited\nHAP', values = c(brewer.pal(9, 'Set1')[9], brewer.pal(9, 'Reds')[7], brewer.pal(9, 'Greens')[7]), drop = FALSE) +

				labs(
					title = paste(opts$pref, paste(chr, NAME, sep = ':'), 'SPRT haplotype inheritance', min.depth),
					subtitle = paste(opts$genmode, proband[[opts$proband]],
						paste0('(#SNP4=', nrow(mm) / 3),
						paste0('na=', format(x = ns[['4a']], digits = 1, nsmall = 1, trim = TRUE)),
						paste0('nb=', format(x = ns[['4b']], digits = 1, nsmall = 1, trim = TRUE)),
						paste0('ff=', format(x = ff, digits = 1, nsmall = 1, trim = TRUE), '%'),
						paste0('pc=', paste(sapply(c(pc_2f, pc_3f$ALL$pc), function(x) format(x = x, digits = 1, nsmall = 1, trim = TRUE)), collapse = '/'), '%'),
						paste0('nc=', paste(sapply(c(nc_2f, pc_3f$ALL$nc), function(x) format(x = x, digits = 1, nsmall = 1, trim = TRUE)), collapse = '/'), '%'),
						paste0('cs=', paste(sapply(c(cs_2f, cs_3f), function(x) format(x = x, digits = 2, nsmall = 2, trim = TRUE)), collapse = '/'), '%'),
						paste0('bs=', format(x = b_score, digits = 2, nsmall = 2, trim = TRUE)),
						paste0('a/b=', paste(opts$alpha, opts$beta, sep = '/'), ')')
					),
					x = '',
					y = ''
				) + theme_minimal() + theme(
					plot.title = element_text(size = rel(1.8), lineheight = 2, vjust = 1, face = 'bold'), plot.subtitle = element_text(size = rel(1.1)),
					legend.title = element_text(size = rel(1.5)), legend.position = 'right', legend.direction = 'vertical',
					axis.title = element_text(size = rel(1.5)), axis.text.x = element_blank(), axis.text.y = element_text(size = rel(1.5))
				)

			# Plotting block arrows
			for(TYPE in c('4a', '4b')){
				x = subset(data.frame(x = 1:nrow(m), subset(m, select = c(SNP.TYPE, DPNI.FWD.HAP, DPNI.REV.HAP))), select = -SNP.TYPE, SNP.TYPE == TYPE)
				y = ifelse(TYPE == '4a', 3.5, 3.4)
				if(nrow(x)){
					mm = sapply(c('FWD', 'REV'), function(X){
						unique(c(.5, x$x[c(which(x[[paste('DPNI', X, 'HAP', sep = '.')]] != 0))], nrow(m) + .5))
					}, simplify = FALSE)
					n = sapply(mm, length, simplify = FALSE)
					m0 = data.frame(
						x0 = mm$FWD[1:(n$FWD-1)],
						x1 = mm$FWD[2:n$FWD], y = rep(y, n$FWD-1)
					)
					m1 = data.frame(
						x1 = mm$REV[1:(n$REV-1)],
						x0 = mm$REV[2:n$REV], y = rep(y-1, n$REV-1)
					)
					p = p +
						geom_segment(data = m0, aes(x = x0, xend = x1, y = y, yend = y), arrow = arrow(length = unit(.5, 'lines'))) +
						geom_segment(data = m1, aes(x = x0, xend = x1, y = y, yend = y), arrow = arrow(length = unit(.5, 'lines'))) +
						geom_segment(data = m0, aes(x = x0, xend = x0, y = y - .04, yend = y + .04)) +
						geom_segment(data = m1, aes(x = x0, xend = x0, y = y - .04, yend = y + .04)) +
						geom_text(data = data.frame(lab = TYPE, y = y), aes(x = 0, y = y, label = lab), hjust = 1, vjust = .5, nudge_x = -1) +
						geom_text(data = data.frame(lab = TYPE, y = y), aes(x = nrow(m), y = y-1, label = lab), hjust = 0, vjust = .5, nudge_x = 1) +
						geom_text(data = data.frame(lab = paste0('(', nrow(m0), ')'), y = y), aes(x = nrow(m), y = y, label = lab), hjust = 0, vjust = .5, nudge_x = 1) +
						geom_text(data = data.frame(lab = paste0('(', nrow(m1), ')'), y = y), aes(x = 0, y = y-1, label = lab), hjust = 1, vjust = .5, nudge_x = -1) +
						geom_segment(data = data.frame(x = 1:nrow(m), y = 1.5, z = as.character(m$SNP.TYPE)),
							aes(x = x, xend = x, y = y - .04, yend = y + .04, colour = factor(z, levels = c('4a', '4b'))), size = 2, alpha = .7) +
						scale_colour_manual(name = 'SNP\nsubtype', values = c(brewer.pal(9, 'Set1')[c(2,4)]), drop = FALSE)
				}
			}

			# Add genomic coordinates
			xt = c(1, c(order(rev(rev(diff(m$POS))[-1]), decreasing = TRUE) + 1)[1:5], nrow(m))
			p = p + geom_text(data = data.frame(x = xt, y = .65, z = m$POS[xt]), aes(x = x, y = y, label = z), size = 2.5, angle = 45, hjust = 1, vjust = .5)

			# Add inter-SNP distance
			p = p + geom_line(data = m, aes(x = 1:nrow(m), y = scales::rescale(x = c(0, diff(POS)), to = c(1-.34, 1+.34))))

			# Pinpoint maternal mutation (if any)
			if(!!length(mut$matmut)){
				if(paste0('chr', mut$matmut$chr) == chr & !is.na(match(NAME, mut$matmut$target))){
					i = c(rev(which(mut$matmut$range[1] >= m$POS))[1], which(mut$matmut$range[2] <= m$POS)[1])
					if(is.na(i[1])){
						i[1] = 1
					}
					if(is.na(i[2])){
						i[2] = nrow(m)
					}
					i = i - .5
					p = p +
						geom_segment(aes(x = i[1], xend = i[2], y = 1.6, yend = 1.6), size = 3) +
						geom_vline(xintercept = i, linetype = 'dashed')
				}
			}
			p = p +
				geom_label(data = nb, aes(x = nrow(m), y = L1), label = paste0('#block=', nb$value), size = 3)
			ggsave(filename = file.path(opts$odir, paste(opts$pref, chr, NAME, 'SNP4_INHERITANCE.png', sep = '_')), plot = p, width = 12, height = 6, device = 'png', dpi = 150)

			m$DPN.HAP = factor(m$DPN.HAP, levels = 0:2, labels = c('?', 'I', 'II'))
			for(a in c('FWD', 'REV')){
				for(b in c('.PROP.RISK.HAP', '.LOWER.QBOUND', '.UPPER.QBOUND')){
					m[[paste0(a, b)]] = format(x = m[[paste0(a, b)]], digits = 3, nsmall = 3, trim = TRUE)
				}
				for(b in c('HAP', 'HAP.SMOOTHED')){
					m[[paste('DPNI', a, b, sep = '.')]] = factor(m[[paste('DPNI', a, b, sep = '.')]], levels = 0:2, labels = c('?', 'I', 'II'))
				}
			}
			m$SNP.TYPE = m$SNP.TYPE.BAK

			write.table(x = subset(m, select = -c(SNP.TYPE.BAK, Fa.ALLELE, Fa.HAP)), file = file.path(opts$odir, paste(opts$pref, chr, NAME, 'SNP4_INHERITANCE.tsv', sep = '_')), sep = "\t", quote = FALSE, row.names = FALSE)

			SPRT.TABLE
		}))
		SPRT.TABLE
	}))
	write.table(x = SPRT.TABLE, file = file.path(opts$odir, paste(opts$pref, 'SPRT.TABLE.tsv', sep = '_')), sep = "\t", quote = FALSE, row.names = FALSE)
}
system(paste('touch', file.path(opts$odir, 'done')))
quit(save = 'no')