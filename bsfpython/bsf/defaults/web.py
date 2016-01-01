"""bsf.defaults.web

A package to centralise web configuration information.
"""

#
# Copyright 2013 - 2016 Michael K. Schuster
#
# Biomedical Sequencing Facility (BSF), part of the genomics core facility
# of the Research Center for Molecular Medicine (CeMM) of the
# Austrian Academy of Sciences and the Medical University of Vienna (MUW).
#
#
# This file is part of BSF Python.
#
# BSF Python is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BSF Python is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with BSF Python.  If not, see <http://www.gnu.org/licenses/>.


import cgi
import datetime
import getpass
import inspect
import urllib


def html_header(title, creator=None, source=None, strict=False):
    """Return the header section of an XHTML 1.0 document.

    @param title: Title element value
    @type title: str
    @param creator: Dublin Core DC.Creator meta field value,
        defaults to the USER environment variable
    @type creator: str
    @param source: Dublin Core DC.Source meta field value,
        defaults to the Python script file path
    @type source: str
    @param strict: XHTML 1.0 Strict or XHTML 1.0 Transitional
    @type strict: bool
    @return: HTML header section as multi-line string
    @rtype: str
    """

    if not creator:
        creator = getpass.getuser()
        # The getpass.getuser method just relies on environment variables,
        # but at least works under Unix and Windows.

    if not source:
        source = inspect.getfile(inspect.currentframe())

    output = str()

    if strict:
        output += '<!DOCTYPE html PUBLIC ' \
                  '"-//W3C//DTD XHTML 1.0 Strict//EN" ' \
                  '"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">\n'
        # output += '<!DOCTYPE HTML PUBLIC ' \
        #           '"-//W3C//DTD HTML 4.01//EN" ' \
        #           '"http://www.w3.org/TR/html4/strict.dtd">\n'
    else:
        output += '<!DOCTYPE html PUBLIC ' \
                  '"-//W3C//DTD XHTML 1.0 Transitional//EN" ' \
                  '"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n'
        # output += '<!DOCTYPE HTML PUBLIC ' \
        #           '"-//W3C//DTD HTML 4.01 Transitional//EN" ' \
        #           '"http://www.w3.org/TR/html4/loose.dtd">\n'

    output += '\n'
    output += '<html xmlns="http://www.w3.org/1999/xhtml">\n'
    output += '<head>\n'
    output += '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1" />\n'
    output += '<link rel="schema.DC" href="http://purl.org/DC/elements/1.0/" />\n'
    output += '<meta name="DC.Creator" content="{}" />\n'.format(creator)
    output += '<meta name="DC.Date" content="{}" />\n'.format(datetime.datetime.now().isoformat())
    output += '<meta name="DC.Source" content="{}" />\n'.format(source)
    output += '<meta name="DC.Title" content="{}" />\n'.format(title)
    output += '<title>{}</title>\n'.format(title)
    output += '</head>\n'
    output += '\n'

    return output


def html_footer():
    """Return the footer section of a XHTML 1.0 document.

    @return: XHTML 1.0 footer section as multi-line string
    @rtype: str
    """

    output = str()

    output += '</html>\n'
    output += '\n'

    return output


# Christoph Bock's ChIP-Seq default track colours.
# All ChIP-Seq factors are indexed in upper case.

chipseq_colours = \
    {
        "OTHER": "0,200,100",
        "HIGH": "57,39,140",
        "MEDIUM": "144,144,144",
        "LOW": "173,0,33",
        "H2": "102,51,204",
        # The following have been prefixed with H3,
        # a distinction needed for shallow peak calling.
        "H3K4ME1": "204,255,51",
        "H3K4ME2": "61,245,0",
        "H3K4ME3": "39,143,68",
        "H3K9ME1": "153,153,153",
        "H3K9ME2": "112,112,112",
        "H3K9ME3": "51,51,51",
        "H3K27ME1": "255,117,71",
        "H3K27ME2": "255,71,10",
        "H3K27ME3": "235,38,39",
        "H3K36ME1": "20,99,255",
        "H3K36ME2": "0,71,214",
        "H3K36ME3": "24,97,174",
        "H3K20ME1": "225,133,35",
        "H3K20ME2": "175,102,24",
        "H3K20ME3": "141,25,28",
        #
        "H4": "204,153,255",
        "AC": "153,0,77",
        #
        "H3K27AC": "153,0,77",
        "H3K56AC": "153,0,77",
        "H4K16AC": "153,0,77",
        #
        "K4ME1": "204,255,51",
        "K4ME2": "61,245,0",
        "K4ME3": "39,143,68",
        "K9ME1": "153,153,153",
        "K9ME2": "112,112,112",
        "K9ME3": "51,51,51",
        "K27ME1": "255,117,71",
        "K27ME2": "255,71,10",
        "K27ME3": "235,38,39",
        "K36ME1": "20,99,255",
        "K36ME2": "0,71,214",
        "K36ME3": "24,97,174",
        "K20ME1": "225,133,35",
        "K20ME2": "175,102,24",
        "K20ME3": "141,25,28",
        "CTCF": "204,153,0",
        "POL": "204,51,77",
        "EZH2": "0,128,153",
        "SUZ12": "128,204,51",
        "RING": "204,204,51",
        "P300": "204,0,0",
        "INPUT": "153,179,204",
        "WCE": "153,179,204"
    }

chipseq_default_factor = 'OTHER'
chipseq_default_colour = '0,0,0'


def get_chipseq_colour(factor=None):
    """Get the web colour for a ChIP-Seq factor.

    If the factor is not in the dictionary, the default colour for factor 'Other' will be returned.
    @param factor: ChIP-Seq factor (e.g. H3K36me3, H3K4me1, ...)
    @type factor: str
    @return: Comma-separated RGB value triplet
    @rtype: str
    """

    if not factor:
        return chipseq_colours[chipseq_default_factor]

    factor_upper = factor.upper()

    if factor_upper in chipseq_colours:
        return chipseq_colours[factor_upper]
    else:
        return chipseq_colours[chipseq_default_factor]


def ucsc_track_url(options_dict, browser_dict=None, track_dict=None, host_name=None):
    """Return a UCSC Genome Browser track URL.

    @param options_dict: Python C{dict} of Python C{str} URL option key value pairs
    @type options_dict: dict
    @param browser_dict: Python C{dict} of Python C{str} browser line key value pairs
    @type browser_dict: dict
    @param track_dict: Python C{dict} of Python C{str} track line (hgct_customText) key value pairs
    @type track_dict: dict
    @param host_name: UCSC Genome Browser host name
    @type host_name: str
    @return: A URL to attach a track to the UCSC Genome Browser
    @rtype: str
    """

    if browser_dict:
        pass

    if track_dict:

        options_dict['hgct_customText'] = 'track'

        keys = track_dict.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            options_dict['hgct_customText'] += ' {}={}'.format(key, track_dict[key])

    if host_name:
        ucsc_name = host_name
    else:
        ucsc_name = 'genome.ucsc.edu'

    primary_url = 'http://{}/cgi-bin/hgTracks?{}'.\
        format(ucsc_name, urllib.urlencode(options_dict))

    return cgi.escape(s=primary_url, quote=True)
