import sys, os, subprocess, shutil
import argparse

# ------------------------------------- #
def main():
    parser = argparse.ArgumentParser(description='change name')
    # parser.add_argument('scriptandargs', help='script and args to launch', nargs='*')
    parser.add_argument('-old', type=str, help='old name', required=True)
    parser.add_argument('-new', type=str, help='new name', required=True)
    parser.add_argument('-dirname', type=str, help='directory to change name', required=True)
    parser.add_argument('--dry', default = False, action="store_true", help='dery run')
    args = vars(parser.parse_args(sys.argv[1:]))
    changename(args)


def is_binary(filename):
    fin = open(filename, 'rb')
    try:
        CHUNKSIZE = 1024
        while 1:
            chunk = fin.read(CHUNKSIZE)
            if '\0' in chunk: # found null byte
                return True
            if len(chunk) < CHUNKSIZE:
                break # done
    finally:
        fin.close()
    return False


class Replace(object):
  def __init__(self, oldname, newname):
    self.old = oldname
    self.new = newname
    self.oldupper = self.old.upper()
    self.oldlower = self.old.lower()
    self.newupper = self.new.upper()
    self.newlower = self.new.lower()
    print ("replace {} ==> {}".format(self.old, self.new))
    print ("replace {} ==> {}".format(self.oldlower, self.oldlower))
    print ("replace {} ==> {}".format(self.oldupper, self.oldupper))
    # sys.exit(1)
  def __call__(self, toto):
    # return  toto.replace(self.old,self.new)
    return  toto.replace(self.oldupper,self.newupper).replace(self.oldlower,self.newlower).replace(self.old,self.new)


def changename(args):
    directory = args['dirname']
    oldname = args['old']
    newname = args['new']
    dry = args['dry']
    exclude = ['.svn', '.DS_Store']
    if not os.path.isdir(directory):
        raise ValueError("directory does not exists " + directory)

    if dry:
        replace = Replace(oldname, newname)
        for root, dirs, files in os.walk(directory, topdown=True):
            dirs[:] = [d for d in dirs if d not in exclude]
            if root.find(".svn") !=-1:
                raise NameError("@@@ ERROR in toto: subfolders must not contain .svn "+ root)
            rootnew = replace(root)
            print('directory {} --> {}'.format(root,rootnew))
            for file in files:
                fileroot = root + '/' + file
                filenew = replace(fileroot)
                print('file {} --> {}'.format(fileroot,filenew))
        return

    backupdirectory = directory + '.old'
    if os.path.isdir(backupdirectory):
        raise ValueError("directory exists " + backupdirectory)
    shutil.copytree(directory, backupdirectory)
    shutil.rmtree(directory)

    replace = Replace(oldname, newname)
    for root, dirs, files in os.walk(backupdirectory, topdown=True):
        dirs[:] = [d for d in dirs if d not in exclude]
        if root.find(".svn") !=-1:
            raise NameError("@@@ ERROR in toto: subfolders must not contain .svn "+ root)
        rootnew = replace(root.replace(backupdirectory,directory))
        print('rootnew', rootnew)
        os.mkdir(rootnew)
        for file in files:
            fileroot = root + '/' + file
            # print 'fileroot: ', fileroot
            # continue
            filenew = replace(fileroot.replace(backupdirectory,directory))
            if fileroot.find('.tgz') !=-1 or  fileroot.find('.png') !=-1 != is_binary(fileroot):
              shutil.copyfile(fileroot, filenew)
              continue
            print('filenew', filenew)
            infile = open(fileroot, 'r')
            outfile = open(filenew, 'w')
            try:
                toto = infile.read()
            except:
                print("cannot read file", fileroot)
            totonew = replace(toto)
            # print 'totonew', totonew
            outfile.write(totonew)
            infile.close()
            outfile.close()

# ------------------------------------- #
if __name__ == '__main__':
    main()
