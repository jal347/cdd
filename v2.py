import pprint
import sys
import os.path
import glob
import time
import logging
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, BeforePosition, AfterPosition
import json

from biothings.utils.common import SubStr, anyfile


class GBFFParser:

    def __init__(self, infile):
        self.infile = infile
        self.in_f = anyfile(self.infile)

    def parse(self):
        """ Uses the dictionary 'gene' created by get_cdd() to create the top layer document for the hub:
        {
            '_id': key of dict
            'regions': value of dict['regions']
            'sites': value of dict['sites']
        }
        see assign_record() for more details

        """

        out_li = []
        gene = self.assign_record()

        for geneid, value in gene.items():
            doc = {
                '_id': geneid,
                'regions': value['regions'],
                'sites': value['sites']
            }
            out_li.append(doc)

        return out_li



    def assign_record(self):
        """ takes a record from the file and creates a document which will be added as a list
        in a dictionary with the geneid as the key.
        The format of each value in this dictionary (see get_site() and get_region() methods for more details):
        {
            'regions': -get_region()
            'sites': -get_site()
        }
        """
        gene = {}

        for rec in SeqIO.parse(self.in_f, 'genbank'):
            gene_id = self.get_geneid(rec)
            region = self.get_region(rec)
            site = self.get_site(rec)
            # checks to see if sites and region contain any values
            if gene_id and (len(site['sites']) > 0 or len(region['regions']) > 0):

                if gene_id in gene.keys():
                    # checks to see if the regions key in get_region(rec) contains at least one value if so append
                    gene[gene_id]['regions'] = gene[gene_id]['regions'] + [region] if len(region['regions']) > 0 else gene[gene_id]['regions']
                    # checks to see if the sties key in get_site(rec) contains at least one value if so append
                    gene[gene_id]['sites'] = gene[gene_id]['sites'] + [site] if len(site['sites']) > 0 else gene[gene_id]['sites']
                else:
                    doc = {
                        'regions': [region] if len(region['regions']) > 0 else [],
                        'sites': [site] if len(site['sites']) > 0 else []
                    }
                    gene[gene_id] = doc
        return gene


    def get_geneid(self, rec):
        """ Gets the gene_id from a record in the file.

        Args: rec: The record is in a GenBank format:
                   https://bioperl.org/formats/sequence_formats/GenBank_sequence_format

        """
        gene_id = None
        gene_feature = [x for x in rec.features if x.type == 'CDS']
        assert len(gene_feature) == 1, '#: {}, id: {}'.format(len(gene_feature), rec.id)
        gene_feature = gene_feature[0]
        db_xref = gene_feature.qualifiers.get('db_xref', None)
        if db_xref:
            x = [x for x in db_xref if x.startswith('GeneID:')]
            if len(x) == 1:
                gene_id = SubStr(x[0], 'GeneID:')
        return gene_id

    def get_region(self, rec):
        """ Creates dict from a record in the file using the regions in the feature section
        {
            'refseq':
            'regions':
                -'id':
                 'location':
                 'name':
                 'note':
        }

        Args: rec: The record is in a GeneBank format:
                   https://bioperl.org/formats/sequence_formats/GenBank_sequence_format
        """
        # make the dict for cdd
        reg = {
            'refseq': rec.id,
            'regions': []
        }
        region_feature = [x for x in rec.features if x.type == 'Region']
        for region in region_feature:

            # get the region location, note, name
            location = region.location
            if location:
                assert len(location.parts) == 1, "There is more than one location, id: {}".format(rec.id)
                # add one because it shows python counting but we want genbank counting
                # https://biopython.org/docs/1.75/api/Bio.SeqFeature.html#Bio.SeqFeature.FeatureLocation.__str__
                location = FeatureLocation(type(location.start)(location.start+1), location.end)
                # has to convert to a string to display fuzzy location "42..>222" if converted to int fuzzy location will
                # not display
                if location.start == location.end:
                    location = str(location.start)
                else:
                    location = str(location.start) + ".." + str(location.end)


            note = region.qualifiers.get('note', None)
            assert note is None or len(note) == 1, "There is more than one cdd note, id:{}".format(rec.id)
            note = note[0] if note else ""

            region_name = region.qualifiers.get('region_name', None)
            assert region_name is None or len(region_name) == 1, "There is more than one region name, id:{}".format(rec.id)
            region_name = region_name[0] if region_name else ""

            # finds the cdd id, checks if there is only one cdd id
            db_xref = region.qualifiers.get('db_xref', None)
            if db_xref:
                x = [x for x in db_xref if x.startswith('CDD:')]
                assert len(x) == 1, "#: {}, id: {}".format(len(x), rec.id)
                cdd_id = SubStr(x[0], 'CDD:')
                doc = {
                    'id': cdd_id,
                    'location': location if location else "",
                    'name': region_name,
                    'note': note
                }
                reg['regions'].append(doc)
            # if no id is found go to next region
            else:
                continue
        return reg

    def get_site(self, rec):
        """ Creates a dictionary from a record in the file for the site section within cdd in the top layer document.
        {
            'refseq':
            'sites':
                - 'id':
                  'location':
                  'type':
                  'note':
        }

        Args: rec: The record is in a GeneBank format:
                   https://bioperl.org/formats/sequence_formats/GenBank_sequence_format
        """

        site = {
            'refseq': rec.id,
            'sites': []
        }

        site_feature = [x for x in rec.features if x.type == 'Site']
        for s in site_feature:

            # get the locations
            sites = ""
            locations = s.location

            if locations:
                locations = s.location.parts
                for location in locations:
                    # add one because it shows python counting but we want genbank counting
                    # https://biopython.org/docs/1.75/api/Bio.SeqFeature.html#Bio.SeqFeature.FeatureLocation.__str__
                    location = FeatureLocation(type(location.start)(location.start + 1), location.end)
                    # has to convert to a string to display fuzzy location "42..>222" if converted to int fuzzy location will
                    # not display
                    if location.start == location.end:
                        location = str(location.start)
                    else:
                        location = str(location.start) + ".." + str(location.end)

                    if sites:
                        sites = sites + "," + location
                    else:
                        sites = location

            # get the site type and notes
            site_type = s.qualifiers.get('site_type', None)
            assert site_type is None or len(site_type) == 1, "There is more than one site_type, id: {}".format(rec.id)
            site_type = site_type[0] if site_type else ""

            note = s.qualifiers.get('note', None)
            assert note is None or len(note) == 1, "There is more than one cdd note, id:{}".format(rec.id)
            note = note[0] if note else ""

            # finds the id, checks if there is only one id
            db_xref = s.qualifiers.get('db_xref', None)
            if db_xref:
                x = [x for x in db_xref if x.startswith('CDD:')]
                assert len(x) == 1, "#: {}, id: {}".format(len(x), rec.id)
                site_id = SubStr(x[0], 'CDD:')
                doc = {
                    'id': site_id,
                    'location': sites,
                    'type': site_type,
                    'note': note
                }
                site['sites'].append(doc)
            # if no id is found go to next site
            else:
                continue
        return site

p = GBFFParser("./protein/zebrafish.1.protein.gpff.gz")

aList = p.parse()

jsonString = json.dumps(aList, indent=2)
jsonFile = open("data/v2fish.json", "w")
jsonFile.write(jsonString)
jsonFile.close()

