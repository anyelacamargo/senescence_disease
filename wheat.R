#save.genome( gdir="./", sdir="./CACHE/", chrs=c("chr1A","chr2A","chr3A","chr4A","chr5A","chr6A","chr7A","chr1B","chr2B","chr3B","chr4B","chr5B","chr6B","chr7B","chr1D","chr2Dx","chr3D","chr4D","chr5D","chr6D","chr7D"), prefix=".wheat", gen=7, haploid=TRUE, mapfile="map.wheat.txt") 

g=load.genome("./CACHE",  chrs=c("chr1A","chr2A","chr3A","chr4A","chr5A","chr6A","chr7A","chr1B","chr2B","chr3B","chr4B","chr5B","chr6B","chr7B","chr1D","chr2Dx","chr3D","chr4D","chr5D","chr6D","chr7D") )

plot.haplotypes <- function( g, subjects=g$subjects, chrs=c("1A","2A","3A","4A","5A","6A","7A","1B","2B","3B","4B","5B","6B","7B","1D","2Dx","3D","4D","5D","6D","7D"), file="haplotypes.pdf", thin=5, subject.file=NULL ) {

    texts = NULL
    subjects.txt = subjects
    if ( !is.null(subject.file)) {
        texts = read.table(subject.file, header=TRUE)
        subjects = texts$CC.Line
        subjects.txt = paste( texts$type, texts$CC.Line )
    }
   
    pdf(file)
    subjects = intersect(subjects,g$subjects)
    s = match( subjects, g$subjects, nomatch=0)
    ns =length(subjects)
    n = max(4,ns) 
    par(mfrow=c(4,1))
    w = g$additive$genome$chromosome %in% chrs
    cat( "nm: ", length(w), "\n")
    w = w[seq(1,length(w),thin)]
    cat( "nm.thin: ", length(w), "\n")
    markers = g$additive$genome$marker[w]
    nm = length(markers)
    bp = g$additive$genome$bp[w]
    bound = which( diff(bp) < 0 )
    print(bound)
    mat = list()
    for ( sub in subjects ) 
        mat[[sub]] = matrix(0,nrow=8,ncol=nm)
    i = 1
    for( m in markers ) {
        d=hdesign(g, m)
        for( j in 1:ns )
            mat[[subjects[j]]][,i] = d[s[j],]
        i = i+1
    }
    
    labels = g$strains
    k = 0
    for ( sub in subjects ) {
        k = k+1

        image( t(mat[[sub]]), axes=FALSE, col=grey(1-(1:10/10)), main=subjects.txt[k])
        
        axis(2,at=seq(0,length(labels)-1)/(length(labels)-1),labels=labels,las=1,cex.axis=0.6)
                                        #		axis(1,at=seq(1,nm)/nm,labels=sprintf("%.2f", g$additive$genome$bp[w]/1.0e6))
        abline( v=bound/nm, col="red")
    }
    
	dev.off()
}	


plot.haplotypes( g,  thin=1, file="wheat.haplotypes.pdf")
