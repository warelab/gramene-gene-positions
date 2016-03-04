'use strict';

var _ = require('lodash');

var levels = ['genome', 'gene', 'transcript', 'protein'];

module.exports = {
  remap: function recursiveRemap(gene, pos, source, dest, transcript_id) {
    var source_idx = levels.indexOf(source);
    var dest_idx = levels.indexOf(dest);
    // invalid level or position < 1 supplied
    if (source_idx === -1 || dest_idx === -1 || pos < 1) {
      return -1;
    }
    if (source === dest) {
      return pos;
    }
    var next_idx = (source_idx < dest_idx) ? source_idx + 1 : source_idx - 1;
    if (next_idx === dest_idx) {
      return calculatePos(gene, pos, source, dest, transcript_id);
    } else { // can't map to dest directly
      var next = levels[next_idx];
      return this.remap(gene, calculatePos(gene, pos, source, next), next, dest, transcript_id);
    }
  }
};

function makeKeys(gene) {
  if (!gene._transcripts) {
    gene._transcripts = _.keyBy(gene.gene_structure.transcripts,'id');
  }
  if (!gene._exons) {
    gene._exons = _.keyBy(gene.gene_structure.exons,'id');
  }
}

function calculatePos(gene, pos, source, dest, transcript_id) {
  switch (source) {
    case 'genome':
      return genomeToGene(gene, pos);
    case 'gene':
      return (dest === 'genome') ? geneToGenome(gene, pos) : geneToTranscript(gene, pos, transcript_id);
    case 'transcript':
      return (dest === 'gene') ? transcriptToGene(gene, pos, transcript_id) : transcriptToProtein(gene, pos, transcript_id);
    case 'protein':
      return proteinToTranscript(gene, pos, transcript_id);
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
  if (pos > (gene.location.end - gene.location.start + 1)) {
    return -1;
  }
  if (gene.location.strand === 1) {
    return gene.location.start + pos - 1;
  }
  return gene.location.end - pos + 1;
}

function geneToTranscript(gene, pos, transcript_id) {
  makeKeys(gene);
  transcript_id = transcript_id || gene.gene_structure.canonical_transcript;
  var exonIds = gene._transcripts[transcript_id].exons;
  var tpos = 0;
  for (var i=0; i<exonIds.length; i++) {
    var exon = gene._exons[exonIds[i]];
    if (pos >= exon.start && pos <= exon.end) {
      return tpos + pos - exon.start + 1;
    }
    tpos += exon.end - exon.start + 1;
  }
  return -1; // gene pos not in the transcript
}

function transcriptToGene(gene, pos, transcript_id) {
  makeKeys(gene);
  transcript_id = transcript_id || gene.gene_structure.canonical_transcript;
  var exonIds = gene._transcripts[transcript_id].exons;
  var pos_in_transcript = 0;
  for (var i=0; i<exonIds.length; i++) {
    var exon = gene._exons[exonIds[i]];
    pos_in_transcript += exon.end - exon.start + 1;
    if (pos <= pos_in_transcript) {
      return exon.end - (pos_in_transcript - pos);
    }
  }
  return -1;
}

function transcriptToProtein(gene, pos, transcript_id) {
  makeKeys(gene);
  transcript_id = transcript_id || gene.gene_structure.canonical_transcript;
  if (!gene._transcripts[transcript_id].hasOwnProperty('cds')) {
    return -1; // no CDS
  }
  var cds = gene._transcripts[transcript_id].cds;
  if (pos > cds.end || pos < cds.start) {
    return -1; // pos not in CDS
  }
  return Math.floor((pos - cds.start) / 3) + 1;
}

function proteinToTranscript(gene, pos, transcript_id) {
  makeKeys(gene);
  transcript_id = transcript_id || gene.gene_structure.canonical_transcript;
  if (!gene._transcripts[transcript_id].hasOwnProperty('cds')) {
    return -1; // no CDS
  }
  var cds = gene._transcripts[transcript_id].cds;
  var tpos = 3 * (pos - 1) + cds.start;
  if (tpos > cds.end || tpos < cds.start) {
    return -1; // position out of range of CDS
  }
  return tpos;
}

