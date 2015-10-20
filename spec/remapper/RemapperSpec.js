describe('Remapper', function () {
  // json response to http://data.gramene.org/genes?idList=AT1G67500,AT4G25880
  // converted into a commonJS module by prepending json doc with
  // `module.exports = `
  var fixture = require('../support/geneFixture.js');
  var remapper = require('../../index.js');

  var geneMinus, genePlus;

  beforeEach(function() {
    geneMinus = fixture.response[0];
    genePlus = fixture.response[1];
  });

  it('should do nothing when mapping between identical levels', function () {
    expect(remapper.remap(geneMinus,25287707,'genome','genome')).toEqual(25287707);
    expect(remapper.remap(geneMinus,10,'gene','gene')).toEqual(10);
    expect(remapper.remap(geneMinus,10,'transcript','transcript')).toEqual(10);
    expect(remapper.remap(geneMinus,10,'protein','protein')).toEqual(10);
    expect(remapper.remap(genePlus,25287707,'genome','genome')).toEqual(25287707);
    expect(remapper.remap(genePlus,10,'gene','gene')).toEqual(10);
    expect(remapper.remap(genePlus,10,'transcript','transcript')).toEqual(10);
    expect(remapper.remap(genePlus,10,'protein','protein')).toEqual(10);
  });
  
  it('should return -1 when given an invalid level', function () {
    expect(remapper.remap(geneMinus,10,'blah','blah')).toEqual(-1);
    expect(remapper.remap(geneMinus,10,'foo','bar')).toEqual(-1);
    expect(remapper.remap(genePlus,10,'blah','blah')).toEqual(-1);
    expect(remapper.remap(genePlus,10,'foo','bar')).toEqual(-1);
  });

  it('should return -1 when given a position < 1', function () {
    expect(remapper.remap(geneMinus,-10,'gene','transcript')).toEqual(-1);
    expect(remapper.remap(genePlus,-10,'gene','transcript')).toEqual(-1);
    expect(remapper.remap(geneMinus,0,'genome','transcript')).toEqual(-1);
    expect(remapper.remap(genePlus,0,'gene','protein')).toEqual(-1);
  });

  it('should return -1 on out-of-range coordinates', function () {
    expect(remapper.remap(geneMinus,100,'genome','gene')).toEqual(-1);
    expect(remapper.remap(geneMinus,999999999,'genome','gene')).toEqual(-1);
    expect(remapper.remap(geneMinus,999999999,'gene','genome')).toEqual(-1);
    expect(remapper.remap(geneMinus,999999999,'gene','transcript')).toEqual(-1);
    expect(remapper.remap(geneMinus,999999999,'transcript','gene')).toEqual(-1);
    expect(remapper.remap(geneMinus,999999999,'transcript','protein')).toEqual(-1);
    expect(remapper.remap(geneMinus,999999999,'protein','transcript')).toEqual(-1);
    expect(remapper.remap(genePlus,100,'genome','gene')).toEqual(-1);
    expect(remapper.remap(genePlus,999999999,'genome','gene')).toEqual(-1);
    expect(remapper.remap(genePlus,999999999,'gene','genome')).toEqual(-1);
    expect(remapper.remap(genePlus,999999999,'gene','transcript')).toEqual(-1);
    expect(remapper.remap(genePlus,999999999,'transcript','gene')).toEqual(-1);
    expect(remapper.remap(genePlus,999999999,'transcript','protein')).toEqual(-1);
    expect(remapper.remap(genePlus,999999999,'protein','transcript')).toEqual(-1);
  });

  it('should return -1 for genomic positions in introns', function() {
    expect(remapper.remap(geneMinus,9100,'gene','transcript')).toEqual(-1);
    expect(remapper.remap(genePlus,500,'gene','transcript')).toEqual(-1);
  });

  it('should remap protein pos 1 to transcript CDS.start', function() {
    expect(remapper.remap(geneMinus,1,'protein','transcript')).toEqual(geneMinus.canonical_transcript.CDS.start);
    expect(remapper.remap(genePlus,1,'protein','transcript')).toEqual(genePlus.canonical_transcript.CDS.start);
  });

  it('should remap protein pos 2 to transcript CDS.start + 3', function() {
    expect(remapper.remap(geneMinus,2,'protein','transcript')).toEqual(geneMinus.canonical_transcript.CDS.start+3);
    expect(remapper.remap(genePlus,2,'protein','transcript')).toEqual(genePlus.canonical_transcript.CDS.start+3);
  });

  it('should remap transcript CDS.start to protein pos 1', function() {
    expect(remapper.remap(geneMinus,geneMinus.canonical_transcript.CDS.start,'transcript','protein')).toEqual(1);
    expect(remapper.remap(genePlus,genePlus.canonical_transcript.CDS.start,'transcript','protein')).toEqual(1);
  });
  
  it('should remap transcript CDS.start+3 to protein pos 2', function() {
    expect(remapper.remap(geneMinus,geneMinus.canonical_transcript.CDS.start+3,'transcript','protein')).toEqual(2);
    expect(remapper.remap(genePlus,genePlus.canonical_transcript.CDS.start+3,'transcript','protein')).toEqual(2);
  });

  it('should remap transcript to gene correctly', function() {
    expect(remapper.remap(geneMinus,1,'transcript','gene')).toEqual(1);
    expect(remapper.remap(genePlus,1,'transcript','gene')).toEqual(1);
    expect(remapper.remap(geneMinus,geneMinus.canonical_transcript.CDS.start,'transcript','gene')).toEqual(224);
  });
  
  it('should remap the protein pos to genomic coord', function() {
    expect(remapper.remap(geneMinus,1,'protein','genome')).toEqual(25296714);
    expect(remapper.remap(genePlus,1,'protein','genome')).toEqual(13155518);
  });
});