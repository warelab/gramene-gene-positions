'use strict';

var levels = ['genome','gene','transcript','protein'];

module.exports = {
  remap: function recursiveRemap(gene, pos, source, dest) {
    var source_idx = levels.indexOf(source);
    var dest_idx = levels.indexOf(dest);
    // invalid level or position < 1 supplied
    if (source_idx === -1 || dest_idx === -1 || pos < 1) {
      return -1;
    }
    if (source === dest) { return pos; }
    var next_idx = (source_idx < dest_idx) ? source_idx + 1: source_idx - 1;
    if (next_idx === dest_idx) {
      return calculatePos(gene, pos, source, dest);
    } else { // can't map to dest directly
      var next = levels[next_idx];
      return this.remap(gene,calculatePos(gene,pos,source,next),next,dest);
    }
  }
};

function calculatePos(gene, pos, source, dest) {
  switch (source) {
    case 'genome':
      return genomeToGene(gene,pos);
    case 'gene':
      return (dest === 'genome') ? geneToGenome(gene,pos) : geneToTranscript(gene,pos);
    case 'transcript':
      return (dest === 'gene') ? transcriptToGene(gene,pos) : transcriptToProtein(gene,pos);
    case 'protein':
      return proteinToTranscript(gene,pos);
  }
}

function genomeToGene(gene, pos) {
  if (pos < gene.location.start || pos > gene.location.end) {
    return -1;
  }
  if (gene.location.strand === 1) {
    return pos - gene.location.start + 1;
  }
  else {
    return gene.location.end - pos + 1;
  }
}

function geneToGenome(gene, pos) {
  if (pos < 1 || pos > gene.location.end - gene.location.start + 1) {
    return -1;
  }
  if (gene.location.strand === 1) {
    return gene.location.start + pos - 1;
  }
  else {
    return gene.location.end - pos + 1;
  }
}

function geneToTranscript(gene, pos) {
  var tpos = 0;
  gene.canonical_transcript.exons.forEach(function(exon) {
    if (pos >= exon.start && pos <= exon.end) {
      tpos += pos - exon.start + 1;
      return tpos;
    }
    tpos += exon.end - exon.start + 1;
  });
  return -1; // gene pos not in the transcript
}

function transcriptToGene(gene, pos) {
  var pos_in_transcript = 0;

  for (var i=0;i<gene.canonical_transcript.exons.length; i++) {
    var exon = gene.canonical_transcript.exons[i];
    pos_in_transcript += exon.end - exon.start + 1;
    if (pos <= pos_in_transcript) {
      return exon.end - (pos_in_transcript - pos);
    }
  }
  return -1;
}

function transcriptToProtein(gene, pos) {
  if (!gene.canonical_transcript.hasOwnProperty('CDS')) {
    return -1; // no CDS
  }
  if (pos > gene.canonical_transcript.CDS.end || pos < gene.canonical_transcript.CDS.start) {
    return -1; // pos not in CDS
  }
  return Math.floor((pos - gene.canonical_transcript.CDS.start)/3) + 1;
}

function proteinToTranscript(gene, pos) {
  if (!gene.canonical_transcript.hasOwnProperty('CDS')) {
    return -1; // no CDS
  }
  var CDS = gene.canonical_transcript.CDS;
  var tpos = 3*(pos-1) + CDS.start;
  if (tpos > CDS.end || tpos < CDS.start) {
    return -1; // position out of range of CDS
  }
  return tpos;
}

