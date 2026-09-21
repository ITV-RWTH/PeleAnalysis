#!/usr/bin/env python3
"""Drop the last stream record from a stream binary directory.

partStream now refuses to produce a stream set that misses a surface node, so
the only way to test how streamBinTubeStats reacts to one is to take a good
directory apart afterwards.  The Header is left untouched: it still claims the
original number of streams, which is what an interrupted or truncated write
looks like from the reader's side.

Usage:  break_stream_bin.py <streamBinDir>
"""

import os
import struct
import sys


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: break_stream_bin.py <streamBinDir>")
    d = sys.argv[1]

    # the Header is text followed by the connectivity as raw bytes, so read it
    # in binary and decode only the leading lines
    header = os.path.join(d, "Header")
    with open(header, "rb") as f:
        f.readline()                       # title
        n_files = int(f.readline())
        f.readline()                       # number of streams
        n_pts = int(f.readline())
        n_comps = int(f.readline())

    # one record: the stream id, then n_comps * n_pts doubles
    record = 4 + 8 * n_comps * n_pts

    # take the record off the last file that holds one
    for i in reversed(range(n_files)):
        name = os.path.join(d, "str_%05d.bin" % i)
        if not os.path.exists(name):
            continue
        with open(name, "rb") as f:
            data = f.read()
        n_in_file = struct.unpack("i", data[:4])[0]
        if n_in_file < 1:
            continue
        with open(name, "wb") as f:
            f.write(struct.pack("i", n_in_file - 1))
            f.write(data[4:len(data) - record])
        print("removed 1 of %d streams from %s" % (n_in_file, name))
        return

    sys.exit("no stream record found in %s" % d)


if __name__ == "__main__":
    main()
