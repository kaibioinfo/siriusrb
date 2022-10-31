require "java"
require "set"
module Sirius
  DEFAULT_PATH = File.expand_path("~/sirius")
  def self.install(path)
    raise "not implemented yet"
  end
end
sirius_home = ENV["SIRIUS_HOME"]
Dir.glob("#{sirius_home}/lib/app/*.jar").each {|jar|
  require(jar)
};nil
require "sirius/depiction.rb"
require "colorize"
module Sirius
  java_import "de.unijena.bioinf.ChemistryBase.ms.Deviation"
end
class Numeric
  def ppm
    Sirius::Deviation.new(self)
  end
  def dalton
    Sirius::Deviation.new(10, self)
  end
end
module Sirius
  class Deviaton
    def to_s
      toString
    end
    def inspect
      toString
    end
  end
  include_package "de.unijena.bioinf.ChemistryBase.chem"
  IAtomContainer = Java::OrgOpenscienceCdkInterfaces::IAtomContainer
  SILENT = Java::OrgOpenscienceCdkSilent::SilentChemObjectBuilder.getInstance
  java_import "de.unijena.bioinf.ChemistryBase.fp.CdkFingerprintVersion"
  java_import "de.unijena.bioinf.MassDecomposer.Chemistry.DecomposerCache"
  DecomposerCacheInstance = DecomposerCache.new(4)

  def decompose(mass, ion: "[M+H]+", dev: 10.ppm, elements: "CHNOPS")
    constr = FormulaConstraints.fromString(elements)
    decomp = DecomposerCacheInstance.get_decomposer(constr.getChemicalAlphabet)
    if ion.nil?
      decomp.decomposeNeutralMassToFormulas(mass, dev, constr)
    else
      ionization = Ion(ion)
      decomp.decomposeToFormulas(mass - ionization.getModification.mass, ionization.ionization, dev, constr)
    end
  end

  def Formula(str)
    return str if str.kind_of?(MolecularFormula)
    MolecularFormula.parseOrThrow(str)
  end
  def show! xs
    if xs.kind_of?(Enumerable)
      opts={}
      ys = xs.map {|x| x._show_molecule(opts)}
    else
      opts={}
      ys = [xs._show_molecule(opts)]
    end
    Utils::Depiction.quick_depict_many(ys, opts)
  end
  def Ion(str)
    PrecursorIonType.fromString(str)
  end
  def SMARTS(query)
    if query.kind_of?(Java::OrgOpenscienceCdkSmarts::SmartsPattern)
      query
    else
      Java::OrgOpenscienceCdkSmarts::SmartsPattern.create(query, SILENT)
    end
  end
  def Mol(obj)
    case obj
    when String
      if obj.downcase.start_with?("inchi=")
        Mol(inchi: obj)
      else
        Mol(smiles: obj)
      end
    when Hash
      inchi = obj[:inchi] || obj["inchi"]
      if inchi
        Java::OrgOpenscienceCdkInchi::InChIGeneratorFactory.getInstance.getInChIToStructure(inchi, SILENT).getAtomContainer
      elsif obj[:sdf] || obj["sdf"] || obj["file"] || obj[:file]
        nm = obj[:sdf] || obj["sdf"]|| obj["file"] || obj[:file]
        reader=nil
        if File.exist?(nm)
          reader = java.io.BufferedReader.new(java.io.FileReader.new("tox21.sdf"))
        else
          reader = java.io.StringReader.new(nm)
        end
        chemreader = Java::OrgOpenscienceCdkIo::ReaderFactory.new().createReader(reader)
        xs=[]
        loop do
          o = Sirius::SILENT.newInstance(Sirius::IAtomContainer)
          break unless chemreader.read(o)
          xs << o
        end
        chemreader.close
        reader.close
        return xs
      else
        smi = obj[:smiles] || obj["smiles"]
        if smi
          Java::OrgOpenscienceCdkSmiles::SmilesParser.new(SILENT).parseSmiles(smi)
        else
          raise "unknown input format #{hash.keys.inspect}"
        end
      end
    when IAtomContainer
      return obj
    end
  end
  IAtomContainer.class_eval do
    def to_s
      Java::OrgOpenscienceCdkSmiles::SmilesGenerator.new(Java::OrgOpenscienceCdkSmiles::SmiFlavor::Generic | Java::OrgOpenscienceCdkSmiles::SmiFlavor::UseAromaticSymbols).create(self)
    end
    def inspect
      to_s
    end
    def ecfp_fingerprint
      return @ecfp_fingerprint if @ecfp_fingerprint
      circ = Java::OrgOpenscienceCdkFingerprint::CircularFingerprinter.new(Java::OrgOpenscienceCdkFingerprint::CircularFingerprinter::CLASS_ECFP6)
      @ecfp_fingerprint = circ.getCountFingerprint(self)
    end
    def pubchem_fingerprint
      return @pubchem_fingerprint if @pubchem_fingerprint
      circ = Java::OrgOpenscienceCdkFingerprint::PubchemFingerprinter.new(SILENT)
      @pubchem_fingerprint = circ.getBitFingerprint(self).asBitSet()
    end
    def maccs_fingerprint
      return @maccs_fingerprint if @maccs_fingerprint
      circ = Java::OrgOpenscienceCdkFingerprint::MACCSFingerprinter.new(SILENT)
      @maccs_fingerprint = circ.getBitFingerprint(self).asBitSet()
    end


    def inchi
      @inchi ||= __calc_inchi__
    end

    def _show_molecule(opts)
      return self
    end

    def __calc_inchi__
      gen = Java::OrgOpenscienceCdkInchi::InChIGeneratorFactory.getInstance.getInChIGenerator(self)
      inchikey = gen.getInchiKey()
      inchi = gen.getInchi()
      InChI.new(inchikey,inchi)
    end
    private :__calc_inchi__

    def tanimoto(other, method=:ecfp)
      other = Mol(other)
      isbinary = false
      a,b = nil
      case method
      when :ecfp, :morgan
        isbinary=false
        a = ecfp_fingerprint()
        b = other.ecfp_fingerprint
      when :pubchem
        isbinary=true
        a = pubchem_fingerprint()
        b = other.pubchem_fingerprint
      when :maccs
        isbinary=true
        a = maccs_fingerprint()
        b = other.maccs_fingerprint
      else
        raise "unknown method #{method}"
      end
      union = 0
      intersection = 0
      if isbinary
        x = a.clone()
        x.and(b)
        intersection = x.cardinality()
        x = a.clone()
        x.or(b)
        union = x.cardinality()
      else
        (0...a.numOfPopulatedbins()).each {|index|
          hashcode = a.getHash(index)
          ca = a.getCount(index)
          cb = b.getCountForHash(hashcode)
          union += [ca, cb].max
          intersection += [ca, cb].min
        }
        (0...b.numOfPopulatedbins()).each {|index|
          hashcode = b.getHash(index)
          ca = a.getCountForHash(hashcode)
          if ca==0
            union += b.getCount(index)
          end
        }
      end
      return intersection.to_f / union.to_f
    end
    def =~ pattern
      if (pattern.kind_of?(String))
        pat = SMARTS(pattern)
        res = pat.match(self)
        if res.length == 0
          return nil
        else
          return Match.new(self, pat, res.to_a)
        end
      elsif pattern.kind_of?(Array)
        pats = pattern.map {|p| SMARTS(p)}.map {|sm| [sm, sm.match(self)]}.select {|sm,r| r.length>0}
        if pats.empty?
          return nil
        else
          return MultiMatch.new(self, pats.map(&:first), pats.map(&:last))
        end
      end
    end 
    def tokenize()
      m = to_s
      atoms = []
      organic = /((?:\()?B|C|N|O|P|S|F|Cl|Br|I|b|c|n|o|p|s)/
      depth = 0
      m.split(/([\[\]])/).each {|token|
        if token == "["
          depth += 1
        elsif token == "]"
          depth -= 1
        elsif depth > 0
          atoms << "[#{token}]" 
        elsif token.empty?
          # ignore
        else
          token.split(organic).select {|x| !x.empty?}.each {|t|
            if t =~ organic
              atoms << t
            else
              atoms.last << t
            end
          }
        end
      }
      atoms
    end
  end
  MolecularFormula.class_eval do
    def to_s
      toString
    end
    def inspect
      toString
    end
    def +(x)
      add(Formula(x))
    end
    def -(x)
      subtract(Formula(x))
    end
    def -@
      negate
    end
    def [](elem)
      numberOf(elem)
    end
  end
  PrecursorIonType.class_eval do
    def to_s
      toString
    end
    def inspect
      toString
    end
  end
  class MultiMatch
    def initialize(molecule, patterns, indizes)
      @molecule = molecule
      @patterns = patterns
      @indizes = indizes
    end

    def show!
      opts={}
      m = [_show_molecule(opts)]
      Utils::Depiction.quick_depict_many(m, opts)
    end
    def _show_molecule(opts)
      done=Set.new
      colors=Sirius::Utils::Depiction::Options2Colors.to_a
      @indizes.zip(colors.cycle).each {|pat,color|
        color=color.first
        chemobjs = Set.new(pat.map {|i| @molecule.getAtom(i)})
        for bond in @molecule.bonds
          if chemobjs.include?(bond.getAtom(0)) && chemobjs.include?(bond.getAtom(1))
            chemobjs.add(bond)
          end
        end

        opts[color] ||= [] 
        opts[color].push(*chemobjs.to_a)
        opts[:glow] = true
      }
      @molecule
    end

    def inspect
      atoms = @molecule.tokenize()
      pal=Sirius::Utils::Depiction::COLOR_NAMES
      indizes={}
      @indizes.each_with_index.map {|i, e|
        i.each {|m| indizes[m] ||= pal[e % pal.size] }
      }
      output=atoms.each_with_index.flat_map {|e,i|
        if indizes[i]
          e.split(/([\(\)])/).map {|x|
            if x=~/[\(\)]/ then x else x.colorize(indizes[i]) end
          }
        else
          e
        end
      }.join()
      output
    end
  end
  class Match
    attr_reader :molecule, :pattern, :indizes
    def initialize(molecule, pattern, indizes)
      @molecule = molecule
      @pattern = pattern
      @indizes = indizes
    end
    def show!
      opts={}
      m = [_show_molecule(opts)]
      Utils::Depiction.quick_depict_many(m, opts)
    end
    def _show_molecule(opts)
      chemobjs = Set.new(indizes.map {|i| @molecule.getAtom(i)})
      for bond in @molecule.bonds
        if chemobjs.include?(bond.getAtom(0)) && chemobjs.include?(bond.getAtom(1))
          chemobjs.add(bond)
        end
      end
      opts[:red] ||= [] 
      opts[:red].push(*chemobjs.to_a)
      opts[:glow] = true
      @molecule
    end
    def inspect
      atoms = @molecule.tokenize()
      indizes = Hash[@indizes.map {|i| [i,i]}]
      output=atoms.each_with_index.flat_map {|e,i|
        if indizes[i]
          e.split(/([\(\)])/).map {|x|
            if x=~/[\(\)]/ then x else x.colorize(:red) end
          }
        else
          e
        end
      }.join()
      def output.inspect
        output
      end
      output
    end

  end
  java_import "de.unijena.bioinf.ChemistryBase.ms.Ms2Experiment"
  java_import "de.unijena.bioinf.babelms.MsIO"
  java_import "de.unijena.bioinf.ChemistryBase.ms.utils.Spectrums"
  java_import "de.unijena.bionf.spectral_alignment.ModifiedCosine"
  java_import "de.unijena.bioinf.sirius.ProcessedInput"
  def MS(file)
    case file
    when String
      MsIO.readExperimentFromFile(java_file(file)).next.mutate() 
    else
      raise "Cannot convert #{file.class} instance to MS spectrum"
    end
  end

  module Ms2Experiment
    # the merged ms2 spectrum
    def preprocessed
      @preprocessed ||= Java::DeUnijenaBioinfSirius::Ms2Preprocessor.new().preprocess(self)
    end
    def cosine(*other, **opts)
      preprocessed.cosine(*other, **opts)
    end
    def modcos(*other, **opts)
      preprocessed.modcos(*other, **opts)
    end
    def revcos(*other, **opts)
      preprocessed.revcos(*other, **opts)
    end
    def ms2
      preprocessed.ms2
    end
    def inspect
      "%s (m/z = %.4f)" % [name, parentmass]
    end
    def parentmass
      ion_mass
    end
  end
  java_import "de.unijena.bionf.spectral_alignment.SpectralSimilarity"
  class SpectralSimilarity
    def shared_peaks
      shardPeaks # I hate this spelling error -_-
    end
    def score
      similarity
    end
    def to_s
      toString
    end
    def inspect
      toString
    end
  end
  class ProcessedInput
    def ms2
      Spectrums.from(getMergedPeaks)
    end
    def parentmass
      getExperimentInformation.ionMass
    end
    def cosine(*other, deviation: Deviation.new(10), sqrt: true, mode: :direct, norm: true, returns: :both)
      case mode
      when :direct, :reversed
        query = Java::DeUnijenaBionfSpectral_alignment::CosineQueryUtils.new(Java::DeUnijenaBionfSpectral_alignment::IntensityWeightedSpectralAlignment.new(deviation))

        queries = other.map {|x| query.createQuery(x.ms2, x.parentmass, false, sqrt)}
        selfquery = query.createQuery(ms2, parentmass, false, sqrt)
        vector = queries.map {|x| mode==:direct ? query.cosineProduct(selfquery, x) : query.cosineProductOfInverse(selfquery, x) }
        case returns
        when :both
          return vector
        when :score, :similarity
          return vector.map(&:score)
        when :shared_peaks
          return vector.map(:shared_peaks)
        else
          raise "unknown return mode #{returns}. Supported are :both, :score, :shared_peaks"
        end
      when :modified
        modcos(*other, deviation: deviation, sqrt: sqrt, mode: mode, norm: norm, returns: returns)
      else
        raise "unknown cosine mode #{mode}, supported are :direct, :reversed, :modified"
      end
    end
    def modcos(*other, deviation: Deviation.new(10), sqrt: true, norm: true, returns: :both)
      otherSpecs = other.map {|x| ModifiedCosine.prepare(x.ms2, x.parentmass, deviation)}
      thisSpec = ModifiedCosine.prepare(ms2, parentmass, deviation)
      vector = other.zip(otherSpecs).map {|x,y|
        ModifiedCosine.new(thisSpec, y, parentmass, x.parentmass, deviation, sqrt ? 0.5 : 1.0).getSimilarity
      }
      if norm
        selfnorm = ModifiedCosine.new(thisSpec, thisSpec, parentmass, parentmass, deviation, sqrt ? 0.5 : 1.0).getScore
        norms = other.zip(otherSpecs).map.map {|x,y| ModifiedCosine.new(y,y, x.parentmass, x.parentmass, deviation, sqrt ? 0.5 : 1.0).getScore}
        vector = vector.zip(norms).map {|v,n| SpectralSimilarity.new(v.score / Math.sqrt(selfnorm*n), v.shared_peaks)}
      end
      case returns
      when :both
        return vector
      when :score, :similarity
        return vector.map(&:score)
      when :shared_peaks
        return vector.map(:shared_peaks)
      else
        raise "unknown return mode #{returns}. Supported are :both, :score, :shared_peaks"
      end
    end
    def revcos(*other, deviation: Deviation.new(10), sqrt: true, norm: true)

    end
    private
    def __cosine__(*other, deviation: Deviation.new(10), sqrt: true, norm: true, returns: :both)
      query = Java::DeUnijenaBionfSpectral_alignment::CosineQueryUtils.new(Java::DeUnijenaBionfSpectral_alignment::IntensityWeightedSpectralAlignment.new(deviation))

      queries = other.map {|x| query.createQuery(x.ms2, x.parentmass, false, sqrt)}
      selfquery = query.createQuery(ms2, parentmass, false, sqrt)



    end
  end

  private
  def java_file(s)
    Java::JavaIo::File.new(s)
  end

end
