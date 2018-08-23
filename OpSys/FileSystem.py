

__all__ = ["FolderStructureTree", "BackupFolder", 
           "getFileModifiedTime", "getFileCreationTime", 
           "serverTimeToLocalTime", "toServerTime"]


import os, time

from datetime import datetime
from dateutil import tz


def getFileModifiedTime(file_path):
    modified_time = datetime.fromtimestamp(os.path.getmtime(file_path))
    return modified_time

def getFileCreationTime(file_path):
    create_time = datetime.fromtimestamp(os.path.getctime(file_path))
    return create_time


def serverTimeToLocalTime(time_server_uct_str):
    zone_server = tz.tzutc()
    zone_local = tz.tzlocal()

    # utc = datetime.utcnow()
    time_server = datetime.strptime(time_server_uct_str, "%Y-%m-%dT%H:%M:%S.%fZ")

    # Tell the datetime object that it's in UTC time zone since 
    # datetime objects are 'naive' by default
    time_server = time_server.replace(tzinfo=zone_server)

    # Convert time zone
    time_local = time_server.astimezone(zone_local)
    time_local = time_local.replace(tzinfo=zone_server)

    return time_local


def toServerTime(time_server_uct_str):
    zone_server = tz.tzutc()
    zone_local = tz.tzlocal()

    # utc = datetime.utcnow()
    time_server = datetime.strptime(time_server_uct_str, "%Y-%m-%dT%H:%M:%S.%fZ")

    # Tell the datetime object that it's in UTC time zone since 
    # datetime objects are 'naive' by default
    time_server = time_server.replace(tzinfo=zone_server)

    return time_server

def dateStrToDatetime(date_str):
    if date_str.find('.')>=0:
        date = datetime.strptime(date_str, "%Y-%m-%d %H:%M:%S.%f")
        return date
    else:
        date = datetime.strptime(date_str, "%Y-%m-%d %H:%M:%S")
        return date


##-----------

def colored(x):
    return('\x1b[31m' + str(x) + '\x1b[0m')



##-------- File structure tree


class FolderStructureTree:
    _IND_ID = 0
    _IND_NAME = 1
    _IND_PARENT = 2
    _IND_CHILDREN = 3
    _IND_FILES = 4
    _IND_METADATA = 5

    _N_ELEMS = 6

    def __init__(self):
        return
        
    def SetDataFile(self, file_path):
        self.dataFileName = file_path

    def SetSyncFolder(self, folder_path):
        self.syncFolderName = folder_path
        
    
    def SetupRootFolderStructure(self, folder_name):
        folder_structure = [None]*self._N_ELEMS
        folder_structure[self._IND_ID] = None
        folder_structure[self._IND_NAME] = folder_name
        folder_structure[self._IND_PARENT] = None
        folder_structure[self._IND_CHILDREN] = []
        folder_structure[self._IND_FILES] = []
        folder_structure[self._IND_METADATA] = None

        #self.SetFolderStructureContents(folder_structure)
        return folder_structure    
        
    def SetFolderStructureContents(self, folder_structure):
        folder = folder_structure[self._IND_NAME]        
                
        folder_content = os.listdir(folder)

        for i in range(len(folder_content)):
            if os.path.isfile(os.path.join(folder, folder_content[i])):
                file_ = [None]*self._N_ELEMS
                file_[self._IND_ID] = None
                file_[self._IND_NAME] = folder_content[i]
                file_[self._IND_PARENT] = folder_structure
                file_[self._IND_CHILDREN] = None
                file_[self._IND_FILES] = None
                file_[self._IND_METADATA] = None
                
                folder_structure[self._IND_FILES].append(file_)
            else:
                #print('---- ', folder_content[i])
                child = [None]*self._N_ELEMS
                child[self._IND_ID] = None
                child[self._IND_NAME] = os.path.join(folder, folder_content[i])
                child[self._IND_PARENT] = folder_structure
                child[self._IND_CHILDREN] = []
                child[self._IND_FILES] = []
                child[self._IND_METADATA] = None
                folder_structure[self._IND_CHILDREN].append(child)
                
        return folder_structure
        
    def SetFolderStructureContents_WalkTree(self, folder_structure):
        folder_list = [folder_structure]
        ind_last = 0
        while ind_last<len(folder_list):
            folder_last = folder_list[ind_last]
            #print(ind_last, folder_last[self._IND_NAME])
            self.SetFolderStructureContents(folder_last)
            sub_folders = folder_last[self._IND_CHILDREN]
            for i in range(len(sub_folders)):
                folder_list.append(sub_folders[i])
                #print(i, sub_folders[i][self._IND_NAME])
            ind_last += 1

    def GetFolderStructure(self, folder_name):
        folder_root = self.SetupRootFolderStructure(folder_name)
        self.SetFolderStructureContents_WalkTree(folder_root)
        return folder_root
    
    def PrintFolderStructure(self, folder_structure, indent=''):
        files = folder_structure[self._IND_FILES]
        print(indent, colored(folder_structure[self._IND_NAME]))
        for i in range(len(files)):
            print(indent, files[i][self._IND_NAME], ' ({})'.format(files[i][self._IND_ID]))
        indent = indent+'   '
        children = folder_structure[self._IND_CHILDREN]
        for i in range(len(children)):
            self.PrintFolderStructure(children[i], indent)
            
    def _GetFolderList_(self, folder_structure, folder_list=[]):
        children = folder_structure[self._IND_CHILDREN]
        for i in range(len(children)):
            folder_list.append(children[i][self._IND_NAME])
        for i in range(len(children)):
            self._GetFolderList_(children[i], folder_list)
        return folder_list

    def GetFolderList(self, folder_structure):
        folder_list = [folder_structure[self._IND_NAME]]
        folder_list = self._GetFolderList_(folder_structure, folder_list)
        return folder_list
        
    def _GetFolderAndIDsList_(self, folder_structure, folder_list=[]):
        children = folder_structure[self._IND_CHILDREN]
        for i in range(len(children)):
            folder_list.append([children[i][self._IND_NAME], children[i][self._IND_ID]])
        for i in range(len(children)):
            self._GetFolderAndIDsList_(children[i], folder_list)
        return folder_list

    def GetFolderAndIDsList(self, folder_structure):
        folder_list = [[folder_structure[self._IND_NAME], folder_structure[self._IND_ID]]]
        folder_list = self._GetFolderAndIDsList_(folder_structure, folder_list)
        return folder_list

    def _GetFileList_(self, folder_structure, file_list=[], folders_ignore=[]):
        folder_name = folder_structure[self._IND_NAME]
        files = folder_structure[self._IND_FILES]
        for i in range(len(files)):
            file_list.append(os.path.join(folder_name, files[i][self._IND_NAME]))
        children = folder_structure[self._IND_CHILDREN]
        for i in range(len(children)):
            if children[i][self._IND_NAME] not in folders_ignore:
                #print(children[i])
                self._GetFileList_(children[i], file_list, folders_ignore)
        return file_list

    def GetFileList(self, folder_structure, folders_ignore=[]):
        file_list = []
        file_list = self._GetFileList_(folder_structure, file_list, folders_ignore=folders_ignore)
        return file_list

    def _GetFileAndIDsList_(self, folder_structure, file_list=[]):
        folder_name = folder_structure[self._IND_NAME]
        files = folder_structure[self._IND_FILES]
        for i in range(len(files)):
            file_name = os.path.join(folder_name, files[i][self._IND_NAME])
            file_id = files[i][self._IND_ID]
            file_list.append([file_name, file_id])
        children = folder_structure[self._IND_CHILDREN]
        for i in range(len(children)):
            self._GetFileAndIDsList_(children[i], file_list)
        return file_list

    def GetFileAndIDList(self, folder_structure):
        file_list = []
        file_list = self._GetFileAndIDsList_(folder_structure, file_list)
        return file_list


    def FindFileByName(self, folder_structure, file_name):
        """ finds the file_name in the folder_structure and return the file structure 
        list ---> file
        file data can be set throught file[self._IND_NAME]=..., file[self._IND_ID]=... 
        """
        root_folder = folder_structure[self._IND_NAME]
        file_name_rel = os.path.relpath(file_name, root_folder)
        folder_parent, file_name_base = os.path.split(file_name_rel)
        assert file_name_base!=''
        folder_list = []
        if folder_parent!='':
            while True:
                folder_parent, folder_name = os.path.split(folder_parent)
                if folder_name!='':
                    folder_list.append(folder_name)
                else:
                    folder_list.reverse()
                    break
        folder_final = folder_structure
        f_ind = 0
        sub_folders_exist = True
        for folder_name in folder_list:
            children = folder_final[self._IND_CHILDREN]
            hit = False
            for child in children:
                if os.path.basename(child[self._IND_NAME])==folder_name:
                    folder_final = child
                    f_ind += 1
                    hit = True
                    break
            if not hit:
                sub_folders_exist = False
                break
        
        if sub_folders_exist:
            files = folder_final[self._IND_FILES]
            ind = None
            for i in range(len(files)):
                if files[i][self._IND_NAME]==file_name_base:
                    ind = i
                    break
            if ind!=None:
                return files[ind]
            else:
                return None
        else:
            return None

                
    def FindFolderByName(self, folder_structure, folder_name):
        """ finds the file_name in the folder_structure and return the file structure 
        list ---> file
        file data can be set throught file[self._IND_NAME]=..., file[self._IND_ID]=... 
        """
        root_folder = folder_structure[self._IND_NAME]
        if os.path.samefile(root_folder, folder_name):
            return folder_structure
        
        folder_name_rel = os.path.relpath(folder_name, root_folder)
        folder_list = []
        folder_parent = folder_name_rel
        if folder_parent!='':
            while True:
                folder_parent, folder_base = os.path.split(folder_parent)
                if folder_base!='':
                    folder_list.append(folder_base)
                else:
                    folder_list.reverse()
                    break
        folder_final = folder_structure
        f_ind = 0
        sub_folders_exist = True
        for folder_base in folder_list:
            children = folder_final[self._IND_CHILDREN]
            hit = False
            for child in children:
                if os.path.basename(child[self._IND_NAME])==folder_base:
                    folder_final = child
                    f_ind += 1
                    hit = True
                    break
            if not hit:
                sub_folders_exist = False
                break
        
        if sub_folders_exist:
            return folder_final
        else:
            return None

    def FindByID(self, folder_structure, _id_):
        """ finds all elements with the given id 
        """
        fs_list = [folder_structure]
        ind_last = 0
        res = []
        while True:
            fs_i = fs_list[ind_last]
            ##skip unwanted folders
            if fs_i[self._IND_ID]==_id_:
                res.append(fs_i)
            fs_files = fs_i[self._IND_FILES]
            for file_i in fs_files:
                if file_i[self._IND_ID]==_id_:
                    res.append(file_i)

            fs_children = fs_i[self._IND_CHILDREN]
            fs_list.extend(fs_children)
            ind_last += 1
            if ind_last>=len(fs_list):
                break
        return res
        


    def SetFileID(self, folder_structure, file_name, file_id):
        """ finds the file_name in the folder_structure and sets its id to file_id
        file_name is a complete file name including the parent folders
        """
        file_node = self.FindFileByName(folder_structure, file_name)
        assert file_node!=None
        file_node[self._IND_ID] = file_id
        return True
        
    def SetFolderID(self, folder_structure, folder_name, folder_id):
        folder_node = self.FindFolderByName(folder_structure, folder_name)
        assert folder_node!=None
        folder_node[self._IND_ID] = folder_id
        return True

    def hasFile(self, folder_structure, file_name):
        file_node = self.FindFileByName(folder_structure, file_name)
        if file_node!=None:
            return True
        else:
            return False

    def hasFolder(self, folder_structure, folder_name):
        folder_node = self.FindFolderByName(folder_structure, folder_name)
        if folder_node!=None:
            return True
        else:
            return False

    def SetFileMetadata(self, folder_structure, file_name, metadata):
        file_node = self.FindFileByName(folder_structure, file_name)
        if file_node!=None:
            file_node[self._IND_METADATA] = metadata
            return True
        else:
            return False
    

    def SetFolderMetadata(self, folder_structure, folder_name, metadata):
        folder_node = self.FindFolderByName(folder_structure, folder_name)
        if folder_node!=None:
            folder_node[self._IND_METADATA] = metadata
            return True
        else:
            return False
    
        
##------- Taking Backup -------------

import shutil
import filecmp

_STR_START = '@$'
_STR_END = '@#&'
_STR_MODIFIED = 'M'
_STR_RENAMED = 'R'
_STR_DELETED = 'D'
_STR_FOLDERSTRUCTURE = 'F'

_DATE_COLON_REP = '='
_DATE_DOT_REP = '#'

class BackupFolder:
    
    def __init__(self, folder_src, folder_dst, folders_ignore=[], file_types_ignore=[]):
        self.folder_src = folder_src
        self.folder_dst = folder_dst
        self.folders_ignore = list(set(self.GetAllSubfolders(folders_ignore)))
        self.fileTypes_ignore = file_types_ignore
        print('Ignoring folders:\n', self.folders_ignore)
        if not os.path.exists(folder_src):
            raise ValueError('Source path does not exist!')
        if not os.path.exists(folder_dst):
            raise ValueError('Destination path does not exist!')
        if not os.path.isdir(folder_src):
            raise ValueError('Source path is not a directory!')
        if not os.path.isdir(folder_dst):
            raise ValueError('Destination path is not a directory!')
            
    def GetAllSubfolders(self, folder_list):
        subs = []
        for folder in folder_list:
            subs += [x[0] for x in os.walk(folder)]
        return subs
    
    def GetFileInfo(self, file_name):
        ## file_name: basename (no folder prefixes)
        ind_start = file_name.find(_STR_START)
        assert ind_start>=0
        ind_start += len(_STR_START)
        ind_end = file_name.find(_STR_END)
        assert ind_end>=0
        info = file_name[ind_start:ind_end]
        name = file_name[ind_end+len(_STR_END):]
        event, date = None, None 
        if info.startswith(_STR_MODIFIED):
            event = 'modified'
            date = info[len(_STR_MODIFIED):]
            date = date.replace(_DATE_COLON_REP, ':').replace(_DATE_DOT_REP, '.')
        elif info.startswith(_STR_DELETED):
            event = 'deleted'
            date = info[len(_STR_DELETED):]
            date = date.replace(_DATE_COLON_REP, ':').replace(_DATE_DOT_REP, '.')
        elif info.startswith(_STR_RENAMED):
            event = 'renamed'
            date = info[len(_STR_RENAMED):]
            date = date.replace(_DATE_COLON_REP, ':').replace(_DATE_DOT_REP, '.')
        elif info.startswith(_STR_FOLDERSTRUCTURE):
            event = 'folderStructure'
            date = info[len(_STR_FOLDERSTRUCTURE):]
            date = date.replace(_DATE_COLON_REP, ':').replace(_DATE_DOT_REP, '.')
        else:
            raise ValueError('Unknown name convention!!')
            
        date = dateStrToDatetime(date)
        return [name, date, event]
        
    def GetFileStamp(self, date, event='modified'):
        """ it returns a unique string containing the date and ... according to 
            a user defined convention
        """
        date = str(date).replace(':', _DATE_COLON_REP).replace('.', _DATE_DOT_REP)
        if event=='modified':
            return _STR_START + _STR_MODIFIED + str(date) + _STR_END
        elif event=='renamed':
            return _STR_START + _STR_RENAMED + str(date) + _STR_END
        elif event=='deleted':
            return _STR_START + _STR_DELETED + str(date) + _STR_END
        elif event=='folderStructure':
            return _STR_START + _STR_FOLDERSTRUCTURE + str(date) + _STR_END
        else:
            raise ValueError()
            
        
    def CopyAll_CloneFolderStruc(self):
        """ it copies all the files and folders in the folder_src to folder_dst
         with the same exact folder structure, the copied files are given names 
         according to the name convention
        """
        src_basename = os.path.basename(self.folder_src)
        src_parent = os.path.split(self.folder_src)[0]
        
        fs_src = FolderStructureTree()
        src_tree = fs_src.GetFolderStructure(self.folder_src)
        
        folder_list = fs_src.GetFolderList(src_tree)
        
        for folder in folder_list:
            dst_name = os.path.relpath(folder, src_parent)
            dst_name = os.path.join(self.folder_dst, dst_name)
            
            if not os.path.exists(dst_name):
                os.makedirs(dst_name)
                
        file_list = fs_src.GetFileList(src_tree)
        for _file_ in file_list:
            
            src_size = os.stat(_file_).st_size
            if src_size>0:
                dst_name = os.path.relpath(_file_, src_parent)
                dst_name = os.path.join(self.folder_dst, dst_name)
                dst_name_folder, dst_name_base = os.path.split(dst_name)
                dst_name_new = self.GetFileStamp(getFileModifiedTime(_file_)) + dst_name_base
                dst_name = os.path.join(dst_name_folder, dst_name_new)
                if not os.path.exists(dst_name):
                    ##check if a file has been renamed by comparing bytes
                    
                    ##otherwise
                    shutil.copyfile(_file_, dst_name)
                    print(_file_, ' --> ', dst_name)
                    dst_size = os.stat(dst_name).st_size
                    print('size: {} ---> {}'.format(src_size, dst_size))
            else:
                print(_file_ + ' : zero size! skipped...')
        return
    
    def CompreFileToFolderContent(self, file_src, folder):
        """ it finds the first identical file inside foder (excluding subfolders)
        """
        folder_content = os.listdir(folder)

        for i in range(len(folder_content)):
            file_dst = os.path.join(folder, folder_content[i])
            if os.path.isfile(file_dst):
                is_identical = filecmp.cmp(file_src, file_dst)
                if is_identical:
                    return file_dst
        return None
    
    
    def CopyAll_Flat(self):
        """ It copies all the file from the source folder and its subfolders 
        to the same folder alongside a file that keeps the source folder structure
        """
        src_basename = os.path.basename(self.folder_src)
        src_parent = os.path.split(self.folder_src)[0]
        
        fs_src = FolderStructureTree()
        src_tree = fs_src.GetFolderStructure(self.folder_src)
        
        folder_list = fs_src.GetFolderList(src_tree)
                        
        file_list = fs_src.GetFileList(src_tree, folders_ignore=self.folders_ignore)
        change = False 
        #print('-'*60)
        #print('file list:\n', file_list)
        
        backup_last = self.GetLatestBackup()
        date_last, fs_tree_last = None, None
        if backup_last!=None:
            date_last, fs_tree_last = backup_last
            print('latest backup date: ', date_last)
            
        for _file_ in file_list:
            src_size = os.stat(_file_).st_size
            
            _file_ext_ = os.path.splitext(_file_)[1]
            if _file_ext_ in self.fileTypes_ignore:
                continue
                
            if src_size>0:
                dst_name = os.path.relpath(_file_, src_parent)
                dst_name = os.path.join(self.folder_dst, dst_name)
                
                dst_name_folder = self.folder_dst
                dst_name_base = os.path.basename(dst_name)
                                
                dst_name_new = self.GetFileStamp(getFileModifiedTime(_file_)) + dst_name_base
                dst_name = os.path.join(dst_name_folder, dst_name_new)
                if not os.path.exists(dst_name):
                    file_clone = self.CompreFileToFolderContent(_file_, self.folder_dst)
                    if file_clone!=None:
                        fs_src.SetFileID(src_tree, _file_, file_clone)
                        # find changes based on modification time -
                        """
                        mtime_file_clone = getFileModifiedTime(file_clone)
                        mtime_file_ = getFileModifiedTime(_file_)
                        if mtime_file_clone<=mtime_file_:
                            os.utime(file_clone)
                            change = True
                            print('!!! NEW COPY: ', _file_)
                        """
                        #check if _file_ is in the latest backup 
                        if fs_tree_last!=None:
                            if not fs_src.hasFile(fs_tree_last, _file_):
                                change = True
                                print(colored('copy or rename: ' + _file_))
                        else:
                            change = True
                            print(colored('copy or rename: ' + _file_))
                    else:
                        shutil.copyfile(_file_, dst_name)
                        print(_file_, ' --> ', dst_name)
                        dst_size = os.stat(dst_name).st_size
                        print('size: {} ---> {}'.format(src_size, dst_size))
                        fs_src.SetFileID(src_tree, _file_, dst_name_new)
                        change = True
                else:
                    #print('file name exists : ', _file_)
                    is_identical = filecmp.cmp(_file_, dst_name)
                    if not is_identical:
                        raise NotImplementedError()
                    
                    fs_src.SetFileID(src_tree, _file_, dst_name_new)
                    if fs_tree_last!=None:
                        if not fs_src.hasFile(fs_tree_last, _file_):
                            change = True
                            print('copy or rename: ', colored(_file_))
                    else:
                        change = True
                        print('copy or rename: ', colored(_file_))

            else:
                print(colored(_file_), ' : zero size! skipped...')

        if change:
            #print('-'*60)       
            #fs_src.PrintFolderStructure(src_tree)
            
            import pickle
            
            fs_tree_name_base = 'FolderStructureTree'
            fs_tree_name = os.path.join(self.folder_dst, fs_tree_name_base)
            file_fs_tree = open(fs_tree_name, 'wb')
            pickle.dump(src_tree, file_fs_tree)
            file_fs_tree.close()
            
            fs_tree_name_base_new = self.GetFileStamp(getFileModifiedTime(fs_tree_name),
                 event='folderStructure') + fs_tree_name_base
            fs_tree_name_new = os.path.join(self.folder_dst, fs_tree_name_base_new)
            os.rename(fs_tree_name, fs_tree_name_new)
        return
        
    
    def FindBackupsSaved(self):
        backups = []
        folder = self.folder_dst
        folder_content = os.listdir(folder)
        
        for i in range(len(folder_content)):
            if os.path.isfile(os.path.join(folder, folder_content[i])):
                fileinfo = self.GetFileInfo(folder_content[i])
                name, date, event = fileinfo
                if event=='folderStructure':
                    name = os.path.join(folder, folder_content[i])
                    import pickle
                    file_ = open(name, 'rb')
                    fs_tree = pickle.load(file_)
                    backups.append([date, fs_tree])
    
        return backups
        
    def FindBackupSorted(self):
        backups = self.FindBackupsSaved()
        dates = [b[0] for b in backups]
        inds = list(range(len(dates)))
        ind_sorted = [i[0] for i in sorted(enumerate(dates), key=lambda x: x[1], reverse=True)]
        backups_sorted = [backups[i] for i in ind_sorted]
        return backups_sorted
        
    def GetLatestBackup(self):
        backups_sorted = self.FindBackupSorted()
        if len(backups_sorted)>0:
            return backups_sorted[0]
        else:
            return None
    
    def RestoreToFolder(self, backup, folder_restore):
        """ backup is a folder tree structure with appropriate file ids
        returned by self.FindBackupsSaved
        """
        #TODO: if the folder name exists inside folder_restore ask for permissions
        ## to override
        folder_src = backup[FolderStructureTree._IND_NAME]
        src_basename = os.path.basename(folder_src)
        src_parent = os.path.split(folder_src)[0]
        
        fs_src = FolderStructureTree()
        src_tree = backup
        
        folder_list = fs_src.GetFolderList(src_tree)
        
        for folder in folder_list:
            dst_name = os.path.relpath(folder, src_parent)
            dst_name = os.path.join(folder_restore, dst_name)
            
            if not os.path.exists(dst_name):
                os.makedirs(dst_name)
                
        file_id_list = fs_src.GetFileAndIDList(src_tree)
        for _file_id_ in file_id_list:
            print(_file_id_)
            _file_, _id_ = _file_id_
            if _id_!=None:
                dst_name = os.path.relpath(_file_, src_parent)
                dst_name = os.path.join(folder_restore, dst_name)
                assert not os.path.exists(dst_name)

                _id_path = os.path.join(self.folder_dst, _id_)
                shutil.copyfile(_id_path, dst_name)
            else:
                print(colored(_file_), ' : ignored...')
        return
    
    
    

