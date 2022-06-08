#!/usr/bin/Rscript
rm(list = ls())
options(warn = -1, bitmapType = 'cairo')
X11.options(type = 'cairo')

rlibs = c('optparse', 'RColorBrewer', 'ggplot2', 'reshape2', 'stringr')
invisible(lapply(rlibs, function(x) suppressMessages(library(x, character.only = TRUE))))

proband = list('AC' = 'Affected child', 'UC' = 'Unaffected child', 'AR' = 'Affected close relative')
gender = list('MM' = 'XY Fetus / XY Proband', 'FF' = 'XX Fetus / XX Proband', 'MF' = 'XY Fetus / XX Proband', 'FM' = 'XX Fetus / XY Proband')

opts = parse_args(object = (args = OptionParser(option_list = list(
	make_option(opt_str = '--pref', type = 'character', help = 'Input files prefix name'),
	make_option(opt_str = '--idir', type = 'character', help = 'Path to input directory'),
	make_option(opt_str = '--odir', type = 'character', help = 'Path to output directory'),
	make_option(opt_str = '--bed', type = 'character', help = 'Path to design BED file'),
	make_option(opt_str = '--gender', type = 'character', help = 'Sexe of fetus and proband (MM, FF, MF or FM)'),
	make_option(opt_str = '--proband', type = 'character', help = 'Proband type (AC, UC or AR)'),
	make_option(opt_str = '--matmut', type = 'character', help = 'Maternal mutation (if any)'),
	make_option(opt_str = '--HTZ', type = 'double', default = .15, help = 'Minimal allele ratio to call an heterozygote site'),
	make_option(opt_str = '--min.SNP', type = 'integer', default = 10, help = 'Minimal number of SNP4 for haplotype inference estimation'),
	make_option(opt_str = '--min.NDP', type = 'integer', default = 8, help = 'Minimal sequencing depth in nuclear dna'),
	make_option(opt_str = '--min.PDP', type = 'integer', default = 15, help = 'Minimal sequencing depth in plasmatic dna'),
	make_option(opt_str = '--min.dist', type = 'integer', default = 50, help = 'Minimal distance between 2 contiguous SNPs for SPRT analysis'),
	make_option(opt_str = '--alpha', type = 'double', default = .001, help = 'alpha risk for SPRT analysis'),
	make_option(opt_str = '--beta', type = 'double', default = .001, help = 'beta risk for SPRT analysis'),
	make_option(opt_str = '--ZFXY', type = 'character', help = 'Path to ZFX/ZFY coverage file in maternal plasmatic DNA')
), usage = 'snptyping.R [options]')))

min.depth = paste0(paste(opts$min.NDP, opts$min.PDP, sep = '/'), 'X')
bedfile = subset(read.table(file = opts$bed, sep = "\t", header = FALSE, as.is = TRUE, col.names = c('CHR', 'START', 'END', 'NAME')), grepl(pattern = 'X', x = CHR))

# Extract mat-/paternal mutation coordinates and target gene(s) (if any)
mut = list()
if(opts$matmut != 'NA'){
	x = unlist(strsplit(x = opts$matmut, split = ':'))
	mut$chr = x[1]
	mut$range = as.numeric(unlist(str_extract_all(string = x[2], pattern = "\\d+")))
	if(length(mut$range) != 2){
		mut$range = rep(mut$range, 2)
	}
	m = subset(bedfile, CHR == mut$chr)
	mut$target = unique(m[c(rev(which(mut$range[1] >= m$START))[1] : which(mut$range[2] <= m$END)[1]), 'NAME'])
}

# PAF/MAF ff estimation
x = read.table(file = file.path(opts$idir, paste0(opts$pref, '_PAF.tsv')), sep = "\t", header = TRUE, as.is = TRUE)
n = paste0(x$Fa, 'RC')
p = 100 * sapply(1:nrow(x), function(i){
	x[[n[i]]][i]
}) / x$TRC
dx = density(p)
xx = dx$x[which.max(dx$y[dx$x < 25])]
ff_paf = 2 * xx
ff_paf.sd = sd(p[p < 100*opts$HTZ])

p = ggplot() +
	geom_density(data = melt(p), aes(value), fill = brewer.pal(9, 'Set1')[2], colour = brewer.pal(9, 'Set1')[2], alpha = .1) +
	geom_vline(xintercept = c(xx - ff_paf.sd, xx, xx + ff_paf.sd), linetype = c('dotted', 'dashed', 'dotted'), alpha = .75) +
	geom_label(data = melt(xx), aes(x = value, y = 0), label = paste0('ff=', round(ff_paf, digits = 1), '%'), size = 5) +
	labs(
		title = paste(opts$pref, 'ff density estimates', min.depth),
		subtitle = paste(gender[[opts$gender]], proband[[opts$proband]], paste0('(#SNP=', length(p), ')')),
		x = 'Paternal Allele Frequency (%)',
		y = 'Density estimates'
	) + theme_minimal() + theme(
		plot.title = element_text(size = rel(1.8), lineheight = 2, vjust = 1, face = 'bold'), plot.subtitle = element_text(size = rel(1.2)),
		legend.title = element_text(size = rel(1.5)), legend.position = 'right', legend.direction = 'vertical',
		axis.title = element_text(size = rel(1.5)), axis.text = element_text(size = rel(1.5))
	)
ggsave(filename = file.path(opts$odir, paste(opts$pref, 'FF_PAF_ESTIMATES.png', sep = '_')), plot = p, width = 10, height = 8, device = 'png', dpi = 150)

x = read.table(file = file.path(opts$idir, paste0(opts$pref, '_MAF.tsv')), sep = "\t", header = TRUE, as.is = TRUE)
MAF = 100 * apply(subset(x, select = c(RRC, ARC)), 1, min) / x$TRC
MAF = MAF[!!MAF]
dx = density(MAF)
xx = dx$x[which.max(dx$y[dx$x < 25])]
ff_maf = 2 * xx
ff_maf.sd = sd(MAF[MAF < 100*opts$HTZ])

p = ggplot() +
	geom_density(data = melt(MAF), aes(value), fill = brewer.pal(9, 'Set1')[2], colour = brewer.pal(9, 'Set1')[2], alpha = .1) +
	geom_vline(xintercept = c(xx - ff_maf.sd, xx, xx + ff_maf.sd), linetype = c('dotted', 'dashed', 'dotted'), alpha = .75) +
	geom_label(data = melt(xx), aes(x = value, y = 0), label = paste0('ff=', round(ff_maf, digits = 1), '%'), size = 5) +
	labs(
		title = paste(opts$pref, 'ff density estimates', min.depth),
		subtitle = paste(gender[[opts$gender]], proband[[opts$proband]], paste0('(#SNP=', length(MAF), ')')),
		x = 'Minor Allele Frequency (%)',
		y = 'Density estimates'
	) + theme_minimal() + theme(
		plot.title = element_text(size = rel(1.8), lineheight = 2, vjust = 1, face = 'bold'), plot.subtitle = element_text(size = rel(1.2)),
		legend.title = element_text(size = rel(1.5)), legend.position = 'right', legend.direction = 'vertical',
		axis.title = element_text(size = rel(1.5)), axis.text = element_text(size = rel(1.5))
	)
ggsave(filename = file.path(opts$odir, paste(opts$pref, 'FF_MAF_ESTIMATES.png', sep = '_')), plot = p, width = 10, height = 8, device = 'png', dpi = 150)

# Add FF estimate from ZFX/ZFY coverage and write to file
covfile = read.table(file = opts$ZFXY, sep = "\t", header = FALSE, as.is = TRUE, col.names = c('CHR', 'START', 'END', 'NAME', 'POS', 'DP'))
ZFXY = list(mu = sapply(split(x = covfile, f = covfile$NAME), function(x){ mean(x$DP) }))
ZFXY$ff = ifelse(test = ZFXY$mu['ZFY'] < 1, yes = NA, no = 200 * ZFXY$mu['ZFY'] / sum(ZFXY$mu))
write.table(x = list(
	PAF = format(x = ff_paf, digits = 1, nsmall = 1, trim = TRUE),
	MAF = format(x = ff_maf, digits = 1, nsmall = 1, trim = TRUE),
	ZFXY = ifelse(test = is.na(ZFXY$ff), yes = NA, no = format(ZFXY$ff, digits = 1, nsmall = 1, trim = TRUE))
	), file = file.path(opts$odir, paste(opts$pref, 'FF.TABLE.tsv', sep = '_')), sep = "\t", quote = FALSE, row.names = FALSE)

# Write ZFY/ZFY and PLASMA (from MAF file) coverage data
write.table(x = lapply(c(ZFXY$mu, PLASMA = mean(x$TRC)), function(x) format(x = x, digits = 1, nsmall = 1, trim = TRUE)),
	file = file.path(opts$odir, paste(opts$pref, 'COV.TABLE.tsv', sep = '_')), sep = "\t", quote = FALSE, row.names = FALSE)

# SPRT analysis
f = ff_maf / 100
q = c(0.5, (1-f)/(2-f), 1/(2-f))
q_bounds = function(dp, q0, q1){
	d = (1-q1) / (1-q0)
	g = (q1 * (1-q0)) / (q0 * (1-q1))
	lapply(list('LOWER.QBOUND' = opts$beta / (1 - opts$alpha), 'UPPER.QBOUND' = (1 - opts$beta) / opts$alpha), function(x){
		(log(x) / dp - log(d)) / log(g)
	})
}

x = read.table(file = file.path(opts$idir, paste0(opts$pref, '_chrX.tsv')), sep = "\t", header = TRUE, as.is = TRUE)
x$BED.TARGET = sapply(x$POS, function(p){
	unique(bedfile$NAME[p >= bedfile$START & p <= bedfile$END])
})

x = split(x = x, f = x$BED.TARGET)

SPRT.TABLE = do.call(rbind, lapply(setdiff(names(x), c('Ethnie', 'ZFX')), function(NAME){
	m = x[[NAME]]
	m = m[c(TRUE, diff(m$POS) >= opts$min.dist), ]
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
			qb = q_bounds(n, q[m$Q0.TYPE[pp[i]]], q[m$Q1.TYPE[pp[i]]])
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
	}
	m = m[order(m$POS), ]

	# Compute %overall identity between FWD and REV DPNI blocks
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

	mm = subset(m, select = c(DPNI.FWD.HAP.SMOOTHED, DPNI.REV.HAP.SMOOTHED))

	names(mm) = c('FWD', 'REV')
	mm$NIPD = 0
	mm$NIPD[i] = mm$FWD[(i = with(mm, FWD == REV))] # Overall DPNI profile ie. FWD/REV cumulative HAP

	mm = data.frame(
		x = rep(1:nrow(mm), 3),
		y = factor(rep(names(mm), each = nrow(mm)), levels = c('NIPD','REV','FWD')),
		z = factor(unlist(mm), levels = c(0:2), labels = c('?', 'I', 'II'))
	)

	nb = melt((n = sapply(c('FWD', 'REV'), function(x){
		length(unique(m[[paste('SPRT', x, 'BLOCK', sep = '.')]]))
	}, simplify = FALSE)))

	# Mean #SNP per block
	ns = mean(sapply(list(m$SPRT.FWD.BLOCK, m$SPRT.REV.BLOCK), function(x) mean(table(x))))

	# Add concordance and block score
	cs_2f = pc_2f
	cs_2f = 1 - logb(100 - cs_2f, base = 100)
	if(length(which(cs_2f > 1)) > 0){
		cs_2f[cs_2f > 1] = 1
	}
	cs_2f[is.na(cs_2f)] = 0
	cs_2f = cs_2f - nc_2f / 100
	cs_2f[which(cs_2f < 0)] = 0

	b_score = min(1 - logb(x =  100 * ns / nrow(m), base = 100), 1)

	# Add SPRT summary to SPRT.TABLE
	SPRT.TABLE = data.frame(
		'GENE' = NAME, # Gene Name
		'#SNP4' = nrow(m), # Total #SNP4
		'#FWD' = n$FWD, # Total #FWD blocks
		'#REV' = n$REV, # Total #REV blocks
		'#MUX' = format(x = ns, digits = 1, nsmall = 1, trim = TRUE), # Mean #SNP4 / block
		'PC' = format(x = pc_2f, digits = 1, nsmall = 1, trim = TRUE),
		'NC' = format(x = nc_2f, digits = 1, nsmall = 1, trim = TRUE),
		'CS' = format(x = cs_2f, digits = 2, nsmall = 2, trim = TRUE),
		'BS' = format(x = b_score, digits = 2, nsmall = 2, trim = TRUE),
	check.names = FALSE)

	p = ggplot(data = mm, aes(x = x, y = y)) +
		geom_tile(aes(fill = z), height = .7) +

		scale_fill_manual(name = 'Inherited\nHAP', values = c(brewer.pal(9, 'Set1')[9], brewer.pal(9, 'Reds')[7], brewer.pal(9, 'Greens')[7]), drop = FALSE) +

		labs(
			title = paste(opts$pref, NAME, 'SPRT haplotype inheritance', min.depth),
			subtitle = paste(gender[[opts$gender]], proband[[opts$proband]],
				paste0('(#SNP4=', nrow(mm) / 3),
				paste0('n=', format(x = ns, digits = 1, nsmall = 1, trim = TRUE)),
				paste0('ff=', format(x = ff_maf, digits = 1, nsmall = 1, trim = TRUE), '%'),
				paste0('pc=', format(x = pc_2f, digits = 1, nsmall = 1, trim = TRUE), '%'),
				paste0('nc=', format(x = nc_2f, digits = 1, nsmall = 1, trim = TRUE), '%'),
				paste0('cs=', format(x = cs_2f, digits = 2, nsmall = 2, trim = TRUE)),
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

	x = data.frame(x = 1:nrow(m), m[, c('DPNI.FWD.HAP', 'DPNI.REV.HAP')])
	y = 3.45
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
		geom_segment(data = m0, aes(x = x0, xend = x0, y = y - .02, yend = y + .02)) +
		geom_segment(data = m1, aes(x = x0, xend = x0, y = y - .02, yend = y + .02))

	# Add genomic coordinates
	xt = c(1, c(order(rev(rev(diff(m$POS))[-1]), decreasing = TRUE) + 1)[1:5], nrow(m))
	p = p + geom_text(data = data.frame(x = xt, y = .65, z = m$POS[xt]), aes(x = x, y = y, label = z), size = 2.5, angle = 45, hjust = 1, vjust = .5)

	# Add inter-SNP distance
	p = p + geom_line(data = m, aes(x = 1:nrow(m), y = scales::rescale(x = c(0, diff(POS)), to = c(1-.34, 1+.34))))

	# Pinpoint maternal mutation (if any)
	if(!!length(mut)){
		if(grepl(pattern = 'X', x = mut$chr) & !is.na(match(NAME, mut$target))){
			i = c(rev(which(mut$range[1] >= m$POS))[1], which(mut$range[2] <= m$POS)[1])
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

	ggsave(filename = file.path(opts$odir, paste(opts$pref, 'chrX', NAME, 'SNP4_INHERITANCE.png', sep = '_')), plot = p, width = 12, height = 6, device = 'png', dpi = 150)

	for(a in c('FWD', 'REV')){
		for(b in c('.PROP.RISK.HAP', '.LOWER.QBOUND', '.UPPER.QBOUND')){
			m[[paste0(a, b)]] = format(x = m[[paste0(a, b)]], digits = 3, nsmall = 3)
		}
		for(b in c('HAP', 'HAP.SMOOTHED')){
			m[[paste('DPNI', a, b, sep = '.')]] = factor(m[[paste('DPNI', a, b, sep = '.')]], levels = 0:2, labels = c('?', 'I', 'II'))
		}
	}
	write.table(x = m, file = file.path(opts$odir, paste(opts$pref, 'chrX', NAME, 'SNP4_INHERITANCE.tsv', sep = '_')), sep = "\t", quote = FALSE, row.names = FALSE)

	return(SPRT.TABLE)
}))
write.table(x = SPRT.TABLE, file = file.path(opts$odir, paste(opts$pref, 'SPRT.TABLE.tsv', sep = '_')), sep = "\t", quote = FALSE, row.names = FALSE)

system(paste('touch', file.path(opts$odir, 'done')))
quit(save = 'no')