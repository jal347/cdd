import sys
import os.path
import glob
import time
import logging
from Bio import SeqIO
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
            'cdd': value of dict
        }
        see get_cdd() for more details

        """

        out_li = []
        gene = self.get_cdd()

        for geneid, value in gene.items():
            doc = {
                '_id': geneid,
                'cdd': value
            }
            out_li.append(doc)

        return out_li

    def get_cdd(self):
        """ takes a record from the file and creates a document which will be added as a list
        in a dictionary with the geneid as the key.
        The format of the document (see get_site() and get_region() methods for more details):
        {
            'refseq': rec.id
            'site': get_site()
            'region': get_region()
        }
        """
        gene = {}

        # goes through each record
        for rec in SeqIO.parse(self.in_f, 'genbank'):

            # calls geneid method use for the key
            gene_id = self.get_geneid(rec)

            # checks to make sure gene_id exists and if a site or region exists
            if gene_id and (len(self.get_site(rec)['site']) > 0 or len(self.get_region(rec)['region']) > 0):
                doc = {
                    'refseq': rec.id,
                    'site': self.get_site(rec)['site'],
                    'region': self.get_region(rec)['region']
                }
                if gene_id in gene.keys():
                    gene[gene_id].append(doc)
                else:
                    gene[gene_id] = [doc]
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
        """ Creates a dictionary from a record in the file for the region section within cdd in the top layer document.
        {
            'region':
                - id:
                  regions:
                  name:
                  note:
        }

        Args: rec: The record is in a GeneBank format:
                   https://bioperl.org/formats/sequence_formats/GenBank_sequence_format
        """
        # make the dict for cdd
        reg = {'region': []}
        region_feature = [x for x in rec.features if x.type == 'Region']
        for region in region_feature:

            # get the region location, note, name
            location = region.location.parts
            assert location is None or len(location) == 1, "There is more than one location, id: {}".format(rec.id)
            if location:
                location = (int(location[0].start + 1), int(location[0].end))

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

                # TODO ASSUMES if id is same notes and region_name is same (NEED TO FIX)
                # this checks if a cdd id is found in the list
                id_index = next((index for (index, d) in enumerate(reg['region']) if d['id'] == cdd_id), None)
                if id_index is not None:
                    if location:
                        # assert note == reg['region'][id_index]['note'], "Notes are not the same, id:{}".format(rec.id)
                        # assert region_name == reg['region'][id_index]['name'], "Region name is not the same, id:{}".format(rec.id)
                        reg['region'][id_index]['regions'].append(location)
                # makes a new entry in CDD if id not found in list
                else:
                    doc = {
                        'id': cdd_id,
                        'regions': [location] if location else [],
                        'name': region_name,
                        'note': note
                    }
                    reg['region'].append(doc)
            # if no id is found go to next region
            else:
                continue
        return reg

    def get_site(self, rec):
        """ Creates a dictionary from a record in the file for the site section within cdd in the top layer document.
        {
            'site':
                - id:
                  sites:
                  type:
        }

        Args: rec: The record is in a GeneBank format:
                   https://bioperl.org/formats/sequence_formats/GenBank_sequence_format
        """

        site = {'site': []}

        site_feature = [x for x in rec.features if x.type == 'Site']
        for s in site_feature:

            # get the locations
            sites = []
            locations = s.location.parts
            if locations:
                for location in locations:
                    location = (int(location.start + 1), int(location.end))
                    sites.append(location)

            site_type = s.qualifiers.get('site_type', None)
            assert site_type is None or len(site_type) == 1, "There is more than one site_type, id: {}".format(rec.id)
            site_type = site_type[0] if site_type else ""

            # finds the id, checks if there is only one id
            db_xref = s.qualifiers.get('db_xref', None)
            if db_xref:
                x = [x for x in db_xref if x.startswith('CDD:')]
                assert len(x) == 1, "#: {}, id: {}".format(len(x), rec.id)
                site_id = SubStr(x[0], 'CDD:')

                # TODO ASSUMES TYPE IS SAME IF ID IS SAME (NEED TO FIX)
                # this checks if a site_id is found in the list
                id_index = next((index for (index, d) in enumerate(site['site']) if d['id'] == site_id), None)
                if id_index is not None:
                    if len(sites) > 0:
                        if site_type != site['site'][id_index]['type']:
                            print("Type is not the same, id:{}".format(rec.id))
                        site['site'][id_index]['sites'] = site['site'][id_index]['sites'] + sites
                # makes a new entry in site if id not found in list
                else:
                    doc = {
                        'id': site_id,
                        'sites': sites,
                        'type': site_type,
                    }
                    site['site'].append(doc)
            # if no id is found go to next site
            else:
                continue
        return site

p = GBFFParser("human.1.protein.gpff.gz")
p.parse()
# aList = p.parse()
#
# jsonString = json.dumps(aList, indent=2)
# jsonFile = open("newdata.json", "w")
# jsonFile.write(jsonString)
# jsonFile.close()
