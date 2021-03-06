## google drive api tools

__all__ = ["GoogleDriveApi", "GDSyncer"]


import httplib2
import os, time, shutil
from send2trash import send2trash

from apiclient import discovery
import oauth2client
from oauth2client import client
from oauth2client import tools

from apiclient import errors
from apiclient import http
from apiclient.http import MediaFileUpload

import argparse
flags = argparse.Namespace(auth_host_name='localhost', auth_host_port=[8080, 8090], 
                           logging_level='ERROR', noauth_local_webserver=False)

from dateutil import tz
#from datetime import tzinfo as tz

from OpSys.FileSystem import FolderStructureTree, BackupFolder, getFileModifiedTime, \
                    getFileCreationTime, serverTimeToLocalTime, toServerTime


class GoogleDriveApi:
    def __init__(self, app_name, client_secret_file, scopes, vbose=False):
        self.app_name = app_name
        self.client_secret_file = client_secret_file
        self.scopes = scopes
        self.vbose = vbose
        
    def get_credentials(self):
        """Gets valid user credentials from storage.

        If nothing has been stored, or if the stored credentials are invalid,
        the OAuth2 flow is completed to obtain the new credentials.

        Returns:
            Credentials, the obtained credential.
        """
        #home_dir = os.path.expanduser('~')
        current_dir = os.getcwd()
        credential_dir = os.path.join(current_dir, 'credentials')
        #credential_dir = './credentials'
        if not os.path.exists(credential_dir):
            os.makedirs(credential_dir)
        credential_path = os.path.join(credential_dir,
                                       (self.app_name+'.json'))
        if self.vbose:
            print('credential_path: ', credential_path)

        store = oauth2client.file.Storage(credential_path)
        credentials = store.get()
        if not credentials or credentials.invalid:
            flow = client.flow_from_clientsecrets(self.client_secret_file, self.scopes)
            flow.user_agent = self.app_name
            if flags:
                credentials = tools.run_flow(flow, store, flags)
            else: # Needed only for compatability with Python 2.6
                credentials = tools.run(flow, store)
            print('Storing credentials to ' + credential_path)
        self.credentials = credentials
        return credentials

    def setupService(self):
        credentials = self.get_credentials()
        self.http = credentials.authorize(httplib2.Http())
        self.service = discovery.build('drive', 'v2', http=self.http)


    def print_file_metadata(self, file_id, pre_str='\t'):
        """Print a file's metadata.

        Args:
            service: Drive API service instance.
            file_id: ID of the file to print metadata for.
        """
        service = self.service
        try:
            file = service.files().get(fileId=file_id).execute()

            print(pre_str+'Title: %s' % file['title'])
            print(pre_str+'MIME type: %s' % file['mimeType'])
            print(pre_str+"createdDate: %s" % file['createdDate'])
            print(pre_str+"modifiedDate: %s" % toServerTime(file['modifiedDate']))
            print(pre_str+"       local: %s" % serverTimeToLocalTime(file['modifiedDate']))
            print(pre_str+"Trashed: %s" % file['labels']['trashed'])
            if file['mimeType']=='application/vnd.google-apps.folder':
                print(pre_str+"type: Folder")
            else:
                print(pre_str+"type: File")
        except errors.HttpError as error:
            print ('An error occurred: %s' % error)

    def get_file_metadata(self, file_id):
        try:
            file_mtd = self.service.files().get(fileId=file_id).execute()
            return file_mtd
        except errors.HttpError as error:
            print ('An error occurred: %s' % error)
            return None

    def get_file_metadata_title(self, file_mtd):
        return file_mtd['title']
    def get_file_metadata_mimeType(self, file_mtd):
        return file_mtd['mimeType']
    def get_file_metadata_createdDate(self, file_mtd):
        return toServerTime(file_mtd['createdDate'])
    def get_file_metadata_modifiedDate(self, file_mtd):
        return toServerTime(file_mtd['modifiedDate'])
    def get_file_metadata_modifiedDate_localTime(self, file_mtd):
        return serverTimeToLocalTime(file_mtd['modifiedDate'])
    def get_file_metadata_isTrashed(self, file_mtd):
        return file_mtd['labels']['trashed']
    def get_file_metadata_isFolder(self, file_mtd):
        if file_mtd['mimeType']=='application/vnd.google-apps.folder':
            return True
        else:
            return False


    def print_file_content(self, file_id):
        """Print a file's content.

        Args:
            service: Drive API service instance.
            file_id: ID of the file.

        Returns:
            File's content if successful, None otherwise.
        """
        service = self.service
        try:
            print(service.files().get_media(fileId=file_id).execute())
        except errors.HttpError as error:
            print('An error occurred: %s' % error)


    def download_file(self, file_id, local_fd):
        """Download a Drive file's content to the local filesystem.

        Args:
            service: Drive API Service instance.
            file_id: ID of the Drive file that will downloaded.
        local_fd: io.Base or file object, the stream that the Drive file's
            contents will be written to.
        """
        service = self.service
        request = service.files().get_media(fileId=file_id)
        media_request = http.MediaIoBaseDownload(local_fd, request)

        print('Download Progress: ', end='')
        while True:
            try:
                download_progress, done = media_request.next_chunk()
            except errors.HttpError as error:
                print('An error occurred: %s' % error)
                return
            if download_progress:
                print(' %d%% ' % int(download_progress.progress() * 100), end='')
            if done:
                print('Download Completed!')
                return

    def insert_file(self, title, description, parent_id, mime_type, filename):
        """Insert new file.

        Args:
            service: Drive API service instance.
            title: Title of the file to insert, including the extension.
            description: Description of the file to insert.
            parent_id: Parent folder's ID.
            mime_type: MIME type of the file to insert.
            filename: Filename of the file to insert.
        Returns:
            Inserted file metadata if successful, None otherwise.
        """
        service = self.service
        media_body = MediaFileUpload(filename, mimetype=mime_type, resumable=True)
        body = {
            'title': title,
            'description': description,
            'mimeType': mime_type
        }
        # Set the parent folder.
        if parent_id:
            body['parents'] = [{'id': parent_id}]

        try:
            request = service.files().insert(body=body, media_body=media_body)
            #file = request.execute()
            ## in chunks
            file = None
            while file==None:
                status, file = request.next_chunk()
                if status:
                    print(" %d%% " % int(status.progress() * 100), end='')
            print()

            return file
        except errors.HttpError as error:
            print('An error occured: %s' % error)
            return None


    def insert_folder(self, title, description='', parent_id=None):
        """Insert new folder.

        Args:
            service: Drive API service instance.
            title: Title of the file to insert, including the extension.
            description: Description of the file to insert.
        Returns:
            Inserted file metadata if successful, None otherwise.
        """
        service = self.service
        body = {
            'title': title,
            'description': description,
            'mimeType': "application/vnd.google-apps.folder"
        }
        # Set the parent folder.
        if parent_id:
            body['parents'] = [{'id': parent_id}]

        try:
            file = service.files().insert(body=body).execute()

            # Uncomment the following line to print the File ID
            # print 'File ID: %s' % file['id']

            return file
        except errors.HttpError as error:
            print('An error occured: %s' % error)
            return None

    def print_children(self, folder_id):
        results = self.service.children().list(folderId=folder_id).execute()
        items = results.get('items', [])
        print('Folder content:')
        children = []
        if not items:
            print('No files found.')
        else:
            for i, item in enumerate(items):
                print('{}: id:({})'.format(i, item['id']))
                self.print_file_metadata(item['id'])


    def get_children(self, folder_id):        
        results = self.service.children().list(folderId=folder_id).execute()
        items = results.get('items', [])
        children = []
        if items:
            for i, item in enumerate(items):
                file_mtd = self.get_file_metadata(item['id'])
                children.append([item['id'], file_mtd])
        return children
        

    def get_subfolders(self, folder_id):        
        results = self.service.children().list(folderId=folder_id).execute()
        items = results.get('items', [])
        subfolds = []
        #print('title:', self.get_file_metadata(folder_id)['title'])
        #print('items:', items)
        if items:
            for i, item in enumerate(items):
                file_mtd = self.get_file_metadata(item['id'])
                if self.get_file_metadata_isFolder(file_mtd):
                    subfolds.append([item['id'], file_mtd])
        return subfolds


    def rename_file(self, file_id, new_title):
        """Rename a file.

        Args:
            service: Drive API service instance.
            file_id: ID of the file to rename.
            new_title: New title for the file.
        Returns:
            Updated file metadata if successful, None otherwise.
        """
        service = self.service
        try:
            file = {'title': new_title}

            # Rename the file.
            updated_file = service.files().patch(
                fileId=file_id,
                body=file,
                fields='title').execute()

            return updated_file
        except errors.HttpError as error:
            print('An error occurred: %s' % error)
            return None


    def update_file(self, file_id, new_title, new_description, new_mime_type,
                    new_filename, new_revision):
        """Update an existing file's metadata and content.

        Args:
            service: Drive API service instance.
            file_id: ID of the file to update.
            new_title: New title for the file.
            new_description: New description for the file.
            new_mime_type: New MIME type for the file.
            new_filename: Filename of the new content to upload.
            new_revision: Whether or not to create a new revision for this file.
        Returns:
            Updated file metadata if successful, None otherwise.
        """
        service = self.service
        try:
            # First retrieve the file from the API.
            file = service.files().get(fileId=file_id).execute()

            # File's new metadata.
            file['title'] = new_title
            file['description'] = new_description
            file['mimeType'] = new_mime_type

            # File's new content.
            media_body = MediaFileUpload(
                new_filename, mimetype=new_mime_type, resumable=True)

            # Send the request to the API.
            update_request = service.files().update(
                fileId=file_id,
                body=file,
                newRevision=new_revision,
                media_body=media_body)
                
            # simple update
            #updated_file = update_request.execute()

            ## update in chunks
            #print('updating file...')
            updated_file = None
            i_ch = 0 #chunk index
            #print('uploading chunk : ', end='')
            while updated_file==None:
                status, updated_file = update_request.next_chunk()
                if status:
                    print(" %d%% " % int(status.progress() * 100), end='')
                #print('{} '.format(i_ch), end='')
                i_ch += 1
            print()

            return updated_file
        except errors.HttpError as error:
            print('An error occurred: %s' % error)
            return None


    def copy_file(self, origin_file_id, copy_title):
        """Copy an existing file.

        Args:
            service: Drive API service instance.
            origin_file_id: ID of the origin file to copy.
            copy_title: Title of the copy.

        Returns:
            The copied file if successful, None otherwise.
        """
        service = self.service
        copied_file = {'title': copy_title}
        try:
            return service.files().copy(
                fileId=origin_file_id, body=copied_file).execute()
        except errors.HttpError as error:
            print('An error occurred: %s' % error)
            return None


    def delete_file(self, file_id):
        """Permanently delete a file, skipping the trash.

        Args:
            service: Drive API service instance.
            file_id: ID of the file to delete.
        """
        service = self.service
        try:
            service.files().delete(fileId=file_id).execute()
        except errors.HttpError as error:
            print('An error occurred: %s' % error)
        

    def retrieve_all_files(self):
        """Retrieve a list of File resources.

        Args:
            service: Drive API service instance.
            Returns:
            List of File resources.
        """
        service = self.service
        result = []
        page_token = None
        while True:
            try:
                param = {}
                if page_token:
                    param['pageToken'] = page_token
                files = service.files().list(**param).execute()

                result.extend(files['items'])
                page_token = files.get('nextPageToken')
                if not page_token:
                    break
            except errors.HttpError as error:
                print('An error occurred: %s' % error)
                break
        return result
        
    def print_files_info(self, n_file_max=100):
        service = self.service
        results = service.files().list(maxResults=n_file_max).execute()
        items = results.get('items', [])
        if not items:
            print('No files found.')
        else:
            print('Files:')
            for i, item in enumerate(items):
                print('{}: {} \n\t id:({}) \n\t type:{}'.format(i, item['title'].encode('utf-8'), 
                    item['id'].encode('utf-8'), item['mimeType'].encode('utf-8')))
        return items
        

    def update_modified_date(self, file_id):
        """Update a file's modified date.

        Args:
            service: Drive API service instance.
            file_id: ID of the file to update the modified date for.

        Returns:
            The updated file if successful, None otherwise.
        """
        service = self.service
        try:
            return service.files().touch(fileId=file_id).execute()
        except errors.HttpError as error:
            print('An error occurred: %s' % error)
        return None
      

    def trash_file(self, file_id):
        """Move a file to the trash.
        
        Args:
            service: Drive API service instance.
            file_id: ID of the file to trash.

        Returns:
            The updated file if successful, None otherwise.
        """
        service = self.service
        try:
            return service.files().trash(fileId=file_id).execute()
        except errors.HttpError as error:
            print('An error occurred: %s' % error)
        return None
      

    def restore_file(self, file_id):
        """Restore a file from the trash.
        
        Args:
            service: Drive API service instance.
            file_id: ID of the file to restore.

        Returns:
            The updated file if successful, None otherwise.
        """
        service = self.service
        try:
            return service.files().untrash(fileId=file_id).execute()
        except errors.HttpError as error:
            print('An error occurred: %s' % error)
        return None
      

    def print_about(self):
        """Print information about the user along with the Drive API settings.
        
        Args:
            service: Drive API service instance.
        """
        service = self.service
        try:
            about = service.about().get().execute()

            print('Current user name: %s' % about['name'])
            print('Root folder ID: %s' % about['rootFolderId'])
            print('Total quota (bytes): %s' % about['quotaBytesTotal'])
            print('Used quota (bytes): %s' % about['quotaBytesUsed'])
            
            return about
        except errors.HttpError as error:
            print('An error occurred: %s' % error)


    def print_change(self, change_id):
        """Print a single Change resource information.

        Args:
            service: Drive API service instance.
            change_id: ID of the Change resource to retrieve.
        """
        service = self.service
        try:
            change = service.changes().get(changeId=change_id).execute()

            print('Changed file ID: %s' % change['fileId'])
            if change['deleted']:
                print('File has been deleted')
            else:
                file = change['file']
                print('Changed file Title: %s' % file['title'])
        except errors.HttpError as error:
            print('An error occurred: %s' % error)


    def retrieve_all_changes(self, start_change_id=None, include_deleted=True, 
            include_subscribed=False, maxResults=None):
        """Retrieve a list of Change resources.

        Args:
            service: Drive API service instance.
            start_change_id: ID of the change to start retrieving subsequent changes
                         from or None.
        Returns:
            List of Change resources.
        """
        service = self.service
        result = []
        page_token = None
        largestChangeId = None
        while True:
            try:
                param = {}
                param['includeDeleted'] = include_deleted
                param['includeSubscribed'] = include_subscribed
                if maxResults!=None:
                    param['maxResults'] = maxResults
                if start_change_id:
                    param['startChangeId'] = start_change_id
                if page_token:
                    param['pageToken'] = page_token
                changes = service.changes().list(**param).execute()
                largestChangeId = changes['largestChangeId']
                print('largestChangeId:', largestChangeId)

                result.extend(changes['items'])
                page_token = changes.get('nextPageToken')
                if not page_token:
                    break
            except errors.HttpError as error:
                print('An error occurred: %s' % error)
                break
        return [largestChangeId, result]

        

    def SetupRootFolderStructure(self, folder_id):
        folder_mtd = self.get_file_metadata(folder_id)
        
        FST = FolderStructureTree
        folder_structure = [None]*FST._N_ELEMS
        folder_structure[FST._IND_ID] = folder_id
        folder_structure[FST._IND_NAME] = self.get_file_metadata_title(folder_mtd)
        folder_structure[FST._IND_PARENT] = None
        folder_structure[FST._IND_CHILDREN] = []
        folder_structure[FST._IND_FILES] = []
        folder_structure[FST._IND_METADATA] = folder_mtd

        return folder_structure    


    def SetFolderStructureContents(self, folder_structure):
        FST = FolderStructureTree
        folder_id = folder_structure[FST._IND_ID]        
                
        id_metadata_list = self.get_children(folder_id)
        for id_metadata in id_metadata_list:
            id, metadata = id_metadata
            if self.get_file_metadata_isFolder(metadata):
                child = [None]*FST._N_ELEMS
                child[FST._IND_ID] = id
                child[FST._IND_NAME] = self.get_file_metadata_title(metadata)
                child[FST._IND_PARENT] = folder_structure
                child[FST._IND_CHILDREN] = []
                child[FST._IND_FILES] = []
                child[FST._IND_METADATA] = metadata
                folder_structure[FST._IND_CHILDREN].append(child)
            else:
                file_ = [None]*FST._N_ELEMS
                file_[FST._IND_ID] = id
                file_[FST._IND_NAME] = self.get_file_metadata_title(metadata)
                file_[FST._IND_PARENT] = folder_structure
                file_[FST._IND_CHILDREN] = None
                file_[FST._IND_FILES] = None
                file_[FST._IND_METADATA] = metadata
                folder_structure[FST._IND_FILES].append(file_)
                
        return folder_structure


    def SetFolderStructureContents_WalkTree(self, folder_structure):
        FST = FolderStructureTree
        folder_list = [folder_structure]
        ind_last = 0
        while ind_last<len(folder_list):
            folder_last = folder_list[ind_last]
            #print(ind_last, folder_last[self._IND_NAME])
            self.SetFolderStructureContents(folder_last)
            sub_folders = folder_last[FST._IND_CHILDREN]
            for i in range(len(sub_folders)):
                folder_list.append(sub_folders[i])
                #print(i, sub_folders[i][self._IND_NAME])
            ind_last += 1

    def GetFolderStructure(self, folder_id):
        folder_root = self.SetupRootFolderStructure(folder_id)
        self.SetFolderStructureContents_WalkTree(folder_root)
        return folder_root

    def PrintFolderStructure(self, folder_structure, indent='  '):
        fs = FolderStructureTree()
        fs.PrintFolderStructure(folder_structure, indent=indent)
        
    def CreateFoldersAndSubfolders(self, folder_local, gd_folder_parent_id='root'):
        FST = FolderStructureTree
        fs = FST().GetFolderStructure(folder_local)
        fs_list = [fs]
        ind_last = 0
        while True:
            fs_i = fs_list[ind_last]
            fold_i_name = os.path.basename(fs_i[FST._IND_NAME])
            fold_i_parent = fs_i[FST._IND_PARENT]
            fold_i_par_id = None
            if fold_i_parent==None:
                fold_i_par_id = gd_folder_parent_id
            else:
                fold_i_par_id = fold_i_parent[FST._IND_ID]
            gdfolder = self.insert_folder(fold_i_name, parent_id=fold_i_par_id)
            fs_i[FST._IND_ID] = gdfolder['id']

            fs_children = fs_i[FST._IND_CHILDREN]
            fs_list.extend(fs_children)
            
            print(fold_i_name, ' created!')

            ind_last += 1
            if ind_last>=len(fs_list):
                break


    def UploadFolder(self, folder_local, gd_folder_parent_id='root', 
            folders_ignore=[], fileTypes_ignore=[], fileEndings_ignore=['~']):
        FST = FolderStructureTree
        fs = FST().GetFolderStructure(folder_local)
        fs_list = [fs]
        ind_last = 0
        while True:
            fs_i = fs_list[ind_last]
            ##skip unwanted folders
            if fs_i[FST._IND_NAME] in folders_ignore:
                print('Folder {} skipped!'.format(fs_i[FST._IND_NAME]))
                ind_last += 1
                if ind_last>=len(fs_list):
                    break
                continue
            fold_i_name = os.path.basename(fs_i[FST._IND_NAME])
            fold_i_parent = fs_i[FST._IND_PARENT]
            fold_i_par_id = None
            if fold_i_parent==None:
                fold_i_par_id = gd_folder_parent_id
            else:
                fold_i_par_id = fold_i_parent[FST._IND_ID]
            gdfolder = self.insert_folder(fold_i_name, parent_id=fold_i_par_id)
            fs_i[FST._IND_ID] = gdfolder['id']
            fs_i[FST._IND_METADATA] = {'modifiedDateServer':toServerTime(gdfolder['modifiedDate'])} 

            print('Folder {} created! id:{}'.format(fold_i_name, gdfolder['id']))
            
            fs_files = fs_i[FST._IND_FILES]
            for file_i in fs_files:
                file_i_title = file_i[FST._IND_NAME]
                ##skip unwanted file types
                _file_i_ext_ = os.path.splitext(file_i_title)[1]
                #print('_file_i_ext_: ', _file_i_ext_)
                if _file_i_ext_ in fileTypes_ignore or file_i_title[-1] in fileEndings_ignore:
                    print('    File {} skipped.'.format(file_i_title))
                    continue
                
                file_i_path = os.path.join(fs_i[FST._IND_NAME], file_i_title)
                ##TODO: set appropriate mime type
                file_i_mime = 'app/my-file'
                gdfile = self.insert_file(title=file_i_title, description='', 
                    parent_id=gdfolder['id'], mime_type=file_i_mime, filename=file_i_path)
                ##print('upload result: ', gdfile)
                if not gdfile:
                    raise ValueError('Error uploading file')
                file_i[FST._IND_ID] = gdfile['id']
                file_i[FST._IND_METADATA] = {'modifiedDateServer':toServerTime(gdfile['modifiedDate']), \
                    'modifiedDateLocal':getFileModifiedTime(file_i_path).replace(tzinfo=tz.tzlocal())} 
                print('    File {} uploaded.'.format(file_i_title))

            fs_children = fs_i[FST._IND_CHILDREN]
            fs_list.extend(fs_children)
            
            ind_last += 1
            if ind_last>=len(fs_list):
                break
            
        return fs
    
    
    def DownloadFolder(self, gd_folder_id, folder_local, overwrite=True):
        """ download gd_folder_id and save it inside folder_local
        """
        if not os.path.exists(folder_local):
            os.makedirs(folder_local)
            
        gd_folder_mtd = self.get_file_metadata(gd_folder_id)
        gdfoldid_list = [[gd_folder_id, gd_folder_mtd, folder_local]]  ##[fold_id, meta data, parent_local] 
        if not self.get_file_metadata_isFolder(gd_folder_mtd):
            raise ValueError()
            
        fst_params_dic = {}
        fold_local_root = None
        
        ind_last = 0
        while True:
            fold_id, fold_mtd, fold_par_local = gdfoldid_list[ind_last]
            children = self.get_children(fold_id)    ##[id, metadata]
            
            ##get fold_id title and create folder local
            fold_title = self.get_file_metadata_title(fold_mtd)
            assert self.get_file_metadata_isFolder(fold_mtd)==True
            fold_local = os.path.join(fold_par_local, fold_title)
            if fold_local_root==None:
                fold_local_root = fold_local
                ##overwrite
                if overwrite and os.path.exists(fold_local_root):
                    resp = input('{} already exists. Overwrite? [y/n]'.format(fold_local_root))
                    if resp=='y':
                        send2trash(fold_local_root)
                        print('{} was deleted!'.format(fold_local_root))
                    else:
                        return None

            if not os.path.exists(fold_local):
                os.makedirs(fold_local)
                
            fst_params_dic[fold_local] = {'type':'folder','id':fold_id, \
                'metadata':{'modifiedDateServer':toServerTime(fold_mtd['modifiedDate'])}}
            
            print('Folder: ', fold_title)
            for f_id, f_mtd in children:
                if self.get_file_metadata_isFolder(f_mtd):
                    ##add to list
                    gdfoldid_list.append([f_id, f_mtd, fold_local])
                    
                else:
                    ##download file
                    f_title = self.get_file_metadata_title(f_mtd)
                    f_local = os.path.join(fold_local, f_title)
                    print('Downloading file :', f_title)
                    self.download_file(f_id, open(f_local, 'wb'))

                    fst_params_dic[f_local] = {'type':'file', 'id':f_id, \
                        'metadata':{'modifiedDateServer':toServerTime(f_mtd['modifiedDate']), \
                        'modifiedDateLocal':getFileModifiedTime(f_local).replace(tzinfo=tz.tzlocal())} }
            ind_last += 1
            if ind_last>=len(gdfoldid_list):
                break
        
        ##create fst (local folder should have been deleted beforehand)
        FST = FolderStructureTree
        fs = FST().GetFolderStructure(fold_local_root)
        fs_list = [fs]
        ind_last = 0
        while True:
            fs_i = fs_list[ind_last]
            f_local = fs_i[FST._IND_NAME]        
            assert f_local in fst_params_dic
            fs_i[FST._IND_ID] = fst_params_dic[f_local]['id']
            fs_i[FST._IND_METADATA] = fst_params_dic[f_local]['metadata']

            fs_files = fs_i[FST._IND_FILES]
            for file_i in fs_files:
                file_i_title = file_i[FST._IND_NAME]                
                file_i_path = os.path.join(fs_i[FST._IND_NAME], file_i_title)
                assert file_i_path in fst_params_dic
                assert fst_params_dic[file_i_path]['type']=='file'
                file_i[FST._IND_ID] = fst_params_dic[file_i_path]['id']
                file_i[FST._IND_METADATA] = fst_params_dic[file_i_path]['metadata']

            fs_children = fs_i[FST._IND_CHILDREN]
            fs_list.extend(fs_children)
            ind_last += 1
            if ind_last>=len(fs_list):
                break
        self.PrintFolderStructure(fs)
        #print('fs:', fs)
        return fs
    

import pickle

class GDSyncer(GoogleDriveApi):
    latestVersion = 1

    def __init__(self, app_name, client_secret_file, scopes, vbose=False):
        """ dataFileName is the file containing the name of the folders to sync
            their respective google drive folder_id, folders and file types to ignore
            etc.
            It can be a list of dictionaraies with each dictionary keeping the 
            information of one (the higher level) synced folder
            It also keeps the time stamp of all the files inside the folder and
            its subfolders
        """
        GoogleDriveApi.__init__(self, app_name, client_secret_file, scopes, vbose)
                
        self.foldersToSync = []
        ##permissions 
        self.askForDownload = True
        self.askToReplaceDL = True
        self.askForUpload = False
        self.askToReplaceUL = False
        
    def SetDataFileName(self, dataFileName):
        self.dataFileName = dataFileName
        
    def SetTempFolder(self, tempFolder):
        self.tempFolder = tempFolder
        
    
    def AddSyncFolder(self, folder_sync, gdfolder_id, folders_ignore=[]
            , fileTypes_ignore=[], fileEndings_ignore=['~']):
        """ Add the folder to datafilename to be synced
            folder: local folder to be synced
            gdfolder_id, gdfolder_name: id and name of the google drive folder 
                to sync to
        """ 
        if gdfolder_id==None:
            gdfolder_id='root'
        self.foldersToSync.append({'foldname':folder_sync, 'gdfold-id':gdfolder_id,\
                'folds-ignore':folders_ignore, 'fileTypes-ignore':fileTypes_ignore, 'fileEndings-ignore':fileEndings_ignore}) 
        return

    

    ##TODO: check if filename exists when downloading (make sure overwriting is ok)
    ##TODO: Track and resume operations, if internet is interrupted
    ##TODO: there might be inconsistency if the file is modified both on server
    ## and locally and it was renamed locally: the server file will always be deleted
    ## and the local file kept
    ##TODO: ignore empty files
    ##TODO: check if 2 files with the same file name are saved in a server folder (validation)
    ##TODO: check if all local and server folders have the same unique file names (validation after failure)
    ##TODO: treat folder renames (server renames are already detected..)
    ##TODO: if a folder id exists on the server but the title has changed (locally.. 
    ## server folder unchanged) we have a local rename. Or if the folder content is 90% identical
    ## to a folder already saved in fst, we have a folder rename.
    ##TODO: if the same new file exists on server and locally, decide which one to keep
    ##TODO: add google drive version in the saved file
    
    ##FIXME: server folder renames are not updated locally
    

    def UpgradeDataFile(self):
        if self.dataFileName!=None and os.path.exists(self.dataFileName):
            file_dataFile = open(self.dataFileName, 'rb')
            dataFileContent = pickle.load(file_dataFile)
            file_dataFile.close()
            dataFileContent_new = None
            if isinstance(dataFileContent, list):
                ##version 0
                ##TODO: remove version 0
                assert False
                dataFileContent_new = {'ver':1, 'data':dataFileContent}
            else:
                ver = dataFileContent['ver']
                if ver<self.latestVersion:
                    raise NotImplementedError()
            if dataFileContent_new!=None:
                file_dataFile = open(self.dataFileName, 'wb')
                pickle.dump(dataFileContent_new, file_dataFile)
                file_dataFile.close()
                print('Data file upgraded to latest version: {}'.format(dataFileContent_new['ver']))
                
            

    def SyncData(self):
        if not os.path.exists(self.dataFileName):
            ##ask what to do :: upload all, download all, or sync based on last 
            ## modification time
            dataFileList = []
            for folderToSyncDic in self.foldersToSync:
                print('Syncing folder:', folderToSyncDic)
                folder_sync = folderToSyncDic['foldname']
                gdfolder_id = folderToSyncDic['gdfold-id']
                folders_ignore = folderToSyncDic['folds-ignore']
                fileTypes_ignore = folderToSyncDic['fileTypes-ignore']
                fileEndings_ignore = folderToSyncDic['fileEndings-ignore']
            
                upload_all = input('Sync file not found. Upload all? [0: Upload all/ \
 1: Download all/ 2:Sync Based on modification time]')
                if upload_all=='0':     ##upload everything
                    fs_tree = self.UploadFolder(folder_local=folder_sync, gd_folder_parent_id=gdfolder_id, 
                        folders_ignore=folders_ignore, fileTypes_ignore=fileTypes_ignore, fileEndings_ignore=fileEndings_ignore)
                    ##save to datafile
                    largestChangeId, result = self.retrieve_all_changes(start_change_id=None, maxResults=1000)
                    dataFileList.append([int(largestChangeId)+1, folderToSyncDic, fs_tree])

                elif upload_all=='1':   ##download all
                    subfolds = self.get_subfolders(gdfolder_id)
                    #print('Choose subfolder to download:', subfolds)
                    print('\n'.join([str(i)+' : '+str(self.get_file_metadata_title(subfolds[i][1])) for i in range(len(subfolds))]))
                    ind = input('Choose subfolder to download: [0-{}]?'.format(len(subfolds)-1))
                    ind = int(ind)
                    assert 0<=ind<len(subfolds)
                    
                    fs_tree = self.DownloadFolder(subfolds[ind][0], folder_sync)
                    if not fs_tree:
                        print('Task interrupted. Choose a different target folder.')
                        return
                    ##save to datafile
                    largestChangeId, result = self.retrieve_all_changes(start_change_id=None, maxResults=1000)
                    dataFileList.append([int(largestChangeId)+1, folderToSyncDic, fs_tree])
                elif upload_all=='2':
                    raise NotImplementedError()
                else:
                    raise ValueError()
                

            self.foldersToSync = []
            dataFileContent = {'ver':self.latestVersion, 'data':dataFileList}
            file_dataFile = open(self.dataFileName, 'wb')
            pickle.dump(dataFileContent, file_dataFile)
            file_dataFile.close()
            print('Done!')
            return True
        else:
            ## open self.dataFileName and check for modifications on server
            ## and on local computer
            self.UpgradeDataFile()
            file_dataFile = open(self.dataFileName, 'rb')
            dataFileContent = pickle.load(file_dataFile)
            file_dataFile.close()
            
            dataFileList = dataFileContent['data']
            
            local_changes_list = self.GetLocalChanges(dataFileList)
            server_changes_list = self.GetServerChanges(dataFileList)
            
            FST = FolderStructureTree
            for i_dfl in range(len(dataFileList)):
                start_change_id, folderToSyncDic, fs_tree = dataFileList[i_dfl]
                local_changes = local_changes_list[i_dfl]
                server_changes = server_changes_list[i_dfl]

                folder_sync = folderToSyncDic['foldname']
                gdfolder_id = folderToSyncDic['gdfold-id']
                folders_ignore = folderToSyncDic['folds-ignore']
                fileTypes_ignore = folderToSyncDic['fileTypes-ignore']
                fileEndings_ignore = folderToSyncDic['fileEndings-ignore']

                for fs_i in local_changes['deleted-folders']:
                    if fs_i not in server_changes['deleted-folders']:
                        fold_i_id = fs_i[FST._IND_ID]
                        fold_i_path = fs_i[FST._IND_NAME]
                        self.trash_file(fold_i_id)
                        print('folder: {} deleted on server.'.format(fold_i_path))
                    if fs_i in fs_i[FST._IND_PARENT][FST._IND_CHILDREN]:
                        fs_i[FST._IND_PARENT][FST._IND_CHILDREN].remove(fs_i)
                
                for fs_i in server_changes['deleted-folders']:
                    if fs_i not in local_changes['deleted-folders']:
                        fold_i_path = fs_i[FST._IND_NAME]
                        if os.path.exists(fold_i_path):
                            send2trash(fold_i_path)
                            print('folder: {} deleted locally.'.format(fold_i_path))
                    if fs_i in fs_i[FST._IND_PARENT][FST._IND_CHILDREN]:
                        fs_i[FST._IND_PARENT][FST._IND_CHILDREN].remove(fs_i)

                for file_i in local_changes['deleted-files']:
                    if file_i not in server_changes['deleted-files']:
                        file_i_id = file_i[FST._IND_ID]
                        self.trash_file(file_i_id)
                        file_i_path = os.path.join(file_i[FST._IND_PARENT][FST._IND_NAME], file_i[FST._IND_NAME])
                        print('file: {} deleted on server.'.format(file_i_path))
                    if file_i in file_i[FST._IND_PARENT][FST._IND_FILES]:
                        file_i[FST._IND_PARENT][FST._IND_FILES].remove(file_i)

                for file_i in server_changes['deleted-files']:
                    if file_i not in local_changes['deleted-files']:
                        file_i_path = os.path.join(file_i[FST._IND_PARENT][FST._IND_NAME], file_i[FST._IND_NAME])
                        if os.path.exists(file_i_path):
                            send2trash(file_i_path)
                            print('file: {} deleted locally.'.format(file_i_path))
                    if file_i in file_i[FST._IND_PARENT][FST._IND_FILES]:
                        file_i[FST._IND_PARENT][FST._IND_FILES].remove(file_i)

                for file_i in local_changes['modified-files']:
                    file_i_id = file_i[FST._IND_ID]
                    file_i_title = file_i[FST._IND_NAME]
                    file_i_path = os.path.join(file_i[FST._IND_PARENT][FST._IND_NAME], file_i_title)

                    file_i_mtd = file_i[FST._IND_METADATA]
                    file_i_mtd_svr = self.get_file_metadata(file_i_id)
                    file_i_ModTime_server_old = file_i_mtd['modifiedDateServer']
                    file_i_ModTime_local_old = file_i_mtd['modifiedDateLocal']
                    
                    file_i_ModTime_local_new = getFileModifiedTime(file_i_path).replace(tzinfo=tz.tzlocal())
                    file_i_ModTime_server_new = toServerTime(file_i_mtd_svr['modifiedDate'])
                    server_changes_modified_files = [server_changes['modified-files'][i][0] for i in range(len(server_changes['modified-files']))] 
                    if file_i not in server_changes_modified_files:
                        updated_gdfile = self.update_file(file_i_id, new_title=file_i_title, new_description='', \
                            new_mime_type=file_i_mtd_svr['mimeType'], new_filename=file_i_path, new_revision=None)
                        file_i[FST._IND_METADATA] = {'modifiedDateServer':toServerTime(updated_gdfile['modifiedDate']), \
                            'modifiedDateLocal':getFileModifiedTime(file_i_path).replace(tzinfo=tz.tzlocal())} 
                        print('modified file: {} uploaded.'.format(file_i_title))
                    else:
                        if file_i_ModTime_local_new>=file_i_ModTime_server_new:
                            updated_gdfile = self.update_file(file_i_id, new_title=file_i_title, new_description='', \
                                new_mime_type=file_i_mtd_svr['mimeType'], new_filename=file_i_path, new_revision=None)
                            file_i[FST._IND_METADATA] = {'modifiedDateServer':toServerTime(updated_gdfile['modifiedDate']), \
                                'modifiedDateLocal':getFileModifiedTime(file_i_path).replace(tzinfo=tz.tzlocal())} 
                            print('modified file: {} uploaded.'.format(file_i_title))
                        else:
                            file_temp_path = os.path.join(self.tempFolder, file_i_title)
                            file_temp = open(file_temp_path, 'wb')
                            self.download_file(file_i_id, local_fd=file_temp)
                            file_temp.close()
                            shutil.copy2(file_temp_path, file_i_path) #copy file and metadata
                            send2trash(file_temp_path)
                            
                            file_i[FST._IND_METADATA] = {'modifiedDateServer':file_i_ModTime_server_new, \
                                'modifiedDateLocal':getFileModifiedTime(file_i_path).replace(tzinfo=tz.tzlocal())} 
                            print('modified file: {} downloaded.'.format(file_i_title))
                    
                for file_i, f_title_new in server_changes['modified-files']:
                    file_i_id = file_i[FST._IND_ID]
                    file_i_title = file_i[FST._IND_NAME]
                    file_i_path = os.path.join(file_i[FST._IND_PARENT][FST._IND_NAME], file_i_title)
                    file_i_path_new = os.path.join(file_i[FST._IND_PARENT][FST._IND_NAME], f_title_new)
                    if not os.path.exists(file_i_path):
                        ##if it was removed and the program was interrupted (to be ckecked for consistency)
                        file_i_path = file_i_path_new

                    file_i_mtd = file_i[FST._IND_METADATA]
                    file_i_mtd_svr = self.get_file_metadata(file_i_id)
                    file_i_ModTime_server_old = file_i_mtd['modifiedDateServer']
                    file_i_ModTime_local_old = file_i_mtd['modifiedDateLocal']
                    
                    file_i_ModTime_local_new = getFileModifiedTime(file_i_path).replace(tzinfo=tz.tzlocal())
                    file_i_ModTime_server_new = toServerTime(file_i_mtd_svr['modifiedDate'])
                    if file_i not in local_changes['modified-files']:
                        file_temp_path = os.path.join(self.tempFolder, file_i_title)
                        file_temp = open(file_temp_path, 'wb')
                        self.download_file(file_i_id, local_fd=file_temp)
                        file_temp.close()
                        shutil.copy2(file_temp_path, file_i_path) #copy file and metadata
                        send2trash(file_temp_path)
                        if file_i_title!=f_title_new:
                            if os.path.exists(file_i_path_new):
                                os.remove(file_i_path_new)
                            os.rename(file_i_path, file_i_path_new)
                            file_i[FST._IND_NAME] = f_title_new

                        file_i[FST._IND_METADATA] = {'modifiedDateServer':file_i_ModTime_server_new, \
                            'modifiedDateLocal':getFileModifiedTime(file_i_path_new).replace(tzinfo=tz.tzlocal())} 
                        if file_i_title!=f_title_new:
                            print('modified file: {} downloaded. (old name: {})'.format(f_title_new, file_i_title))
                        else:
                            print('modified file: {} downloaded.'.format(file_i_title))
                            
                            
                for fold_i, f_title_new, f_modDate in server_changes['modified-folders']:
                    fold_i_id = fold_i[FST._IND_ID]
                    fold_i_path = fold_i[FST._IND_NAME]
                    ##TODO:ensure folder names ending by an slash do not create problems
                    fold_i_parent = os.path.dirname(os.path.normpath(fold_i_path))
                    fold_i_path_new = os.path.join(fold_i_parent, f_title_new)
                    if os.path.normpath(fold_i_path)!=os.path.normpath(fold_i_path_new):
                        ##Folder rename on server
                        print('Folder rename: {} --> {}'.format(fold_i_path, fold_i_path_new))
                        if not os.path.exists(fold_i_path_new):
                            os.rename(fold_i_path, fold_i_path_new)
                            print('Folder successfully renamed.')
                            fold_i[FST._IND_METADATA] = {'modifiedDateServer':toServerTime(f_modDate)}
                            fold_i[FST._IND_NAME] = fold_i_path_new
                        else:
                            print('Error: {} already exists..'.format(fold_i_path_new))


                for fs_i in local_changes['new-folders']:
                    fold_local = fs_i[FST._IND_NAME]
                    fs_parent = FST().FindFolderByName(fs_tree, fs_i[FST._IND_PARENT][FST._IND_NAME])
                    assert fs_parent!=None
                    fs_uploaded = self.UploadFolder(fold_local, gd_folder_parent_id=fs_parent[FST._IND_ID], 
                        folders_ignore=folders_ignore, fileTypes_ignore=fileTypes_ignore, fileEndings_ignore=fileEndings_ignore)
                    fs_uploaded[FST._IND_PARENT] = fs_parent
                    fs_parent[FST._IND_CHILDREN].append(fs_uploaded)
                    print('new folder: {} created on server.'.format(fold_local))

                for file_i in local_changes['new-files']:
                    ##TODO: ckeck if the new file already exists on server
                    file_i_id = file_i[FST._IND_ID]
                    file_i_title = file_i[FST._IND_NAME]
                    file_i_path = os.path.join(file_i[FST._IND_PARENT][FST._IND_NAME], file_i_title)
                    print('new file: {}'.format(file_i_path))
                    file_i_parent = FST().FindFolderByName(fs_tree, file_i[FST._IND_PARENT][FST._IND_NAME])
                    assert file_i_parent!=None

                    file_i_mime = 'app/my-file'
                    gdfile = self.insert_file(title=file_i_title, description='', 
                        parent_id=file_i_parent[FST._IND_ID], mime_type=file_i_mime, filename=file_i_path)
                    if not gdfile:
                        raise ValueError('Error uploading file')
                    file_i[FST._IND_ID] = gdfile['id']
                    file_i[FST._IND_METADATA] = {'modifiedDateServer':toServerTime(gdfile['modifiedDate']), \
                        'modifiedDateLocal':getFileModifiedTime(file_i_path).replace(tzinfo=tz.tzlocal())} 
                    file_i[FST._IND_PARENT] = file_i_parent
                    file_i_parent[FST._IND_FILES].append(file_i)
                    print('new file {} uploaded. '.format(file_i_title))


                for f_id, f_title, f_modDate, p_local in server_changes['new-folders']:
                    f_localpath = os.path.join(p_local[FST._IND_NAME], f_title)
                    fst_new = self.DownloadFolder(gd_folder_id=f_id, folder_local=p_local[FST._IND_NAME], overwrite=True)
                    p_local[FST._IND_CHILDREN].append(fst_new)
                    fst_new[FST._IND_PARENT] = p_local
                    fst_new[FST._IND_METADATA] = {'modifiedDateServer':toServerTime(f_modDate)}
                    print('new folder {} downloaded. '.format(f_localpath))

                for f_id, f_title, f_modDate, p_local  in server_changes['new-files']:
                    file_temp_path = os.path.join(self.tempFolder, f_title)
                    file_temp = open(file_temp_path, 'wb')
                    self.download_file(f_id, local_fd=file_temp)
                    file_temp.close()

                    f_localpath = os.path.join(p_local[FST._IND_NAME], f_title)
                    if os.path.exists(f_localpath):
                        ##TODO: check which one is more recent
                        overwrite = input('New file {} already exists. Overwrite?[y/n]'.format(f_localpath))
                        if overwrite=='y':
                            shutil.copy2(file_temp_path, f_localpath) #copy file and metadata
                            send2trash(file_temp_path)
                            print('new file {} downloaded. '.format(f_localpath))
                            f_local = FST().FindFileByName(p_local, f_localpath)
                            assert f_local!=None
                            f_local[FST._IND_ID] = f_id
                            f_local[FST._IND_METADATA] = {'modifiedDateServer':toServerTime(f_modDate), \
                            'modifiedDateLocal':getFileModifiedTime(f_localpath).replace(tzinfo=tz.tzlocal())} 
                        elif overwrite=='n':
                            ##delete the server version
                            self.trash_file(f_id)
                            print('file: {} deleted on server.'.format(f_title))
                            
                            ##FIXME: local file is automatically uploaded in the local_changes[new-files] wether y/n
                            
                        
                            
                            
                    else:
                        shutil.copy2(file_temp_path, f_localpath) #copy file and metadata
                        send2trash(file_temp_path)
                        print('new file {} downloaded. '.format(f_localpath))
                        f_local = [None]*FST._N_ELEMS
                        f_local[FST._IND_ID] = f_id
                        f_local[FST._IND_NAME] = f_title
                        f_local[FST._IND_PARENT] = p_local
                        f_local[FST._IND_METADATA] = {'modifiedDateServer':toServerTime(f_modDate), \
                        'modifiedDateLocal':getFileModifiedTime(f_localpath).replace(tzinfo=tz.tzlocal())} 
                        p_local[FST._IND_FILES].append(f_local)
                    
                    
                
                largestChangeId, result = self.retrieve_all_changes(dataFileList[i_dfl][0]+10000, maxResults=1000)    
                dataFileList[i_dfl][0] = int(largestChangeId)+1

            dataFileContent = {'ver':self.latestVersion, 'data':dataFileList}
            file_dataFile = open(self.dataFileName, 'wb')
            pickle.dump(dataFileContent, file_dataFile)
            file_dataFile.close()
            print('Done!')
            return True

                        
    def GetLocalChanges(self, dataFileList):
        FST = FolderStructureTree
        local_changes_list = [None]*len(dataFileList)
        for i in range(len(dataFileList)):
            local_changes_list[i] = {'deleted-folders':[], 'deleted-files':[], 
            'modified-folders':[], 'modified-files':[], 
            'new-folders':[], 'new-files':[],
            'renamed-folders':[], 'renamed-files':[],
            'ignored-folders':[], 'ignored-files':[]}

        for i_dfl in range(len(dataFileList)):
            start_change_id, folderToSyncDic, fs_tree = dataFileList[i_dfl]
            local_changes = local_changes_list[i_dfl]
            
            folder_sync = folderToSyncDic['foldname']
            gdfolder_id = folderToSyncDic['gdfold-id']
            folders_ignore = folderToSyncDic['folds-ignore']
            fileTypes_ignore = folderToSyncDic['fileTypes-ignore']
            fileEndings_ignore = folderToSyncDic['fileEndings-ignore']
    
            fs_list = [fs_tree]
            ind_last = 0
            while True:
                fs_i = fs_list[ind_last]
                ##skip unwanted folders
                if fs_i[FST._IND_NAME] in folders_ignore:
                    local_changes['ignored-folders'].append(fs_i)
                    ind_last += 1
                    if ind_last>=len(fs_list):
                        break
                    continue
                
                fold_i_path = fs_i[FST._IND_NAME]
                fold_i_name = os.path.basename(fold_i_path)
                fold_i_id = fs_i[FST._IND_ID]
                
                if not os.path.exists(fold_i_path):
                    local_changes['deleted-folders'].append(fs_i)

                    ind_last += 1
                    if ind_last>=len(fs_list):
                        break
                    continue

                fs_files = fs_i[FST._IND_FILES]
                for file_i in fs_files:
                    file_i_title = file_i[FST._IND_NAME]
                    ##skip unwanted file types
                    _file_i_ext_ = os.path.splitext(file_i_title)[1]
                    if _file_i_ext_ in fileTypes_ignore  or file_i_title[-1] in fileEndings_ignore:
                        local_changes['ignored-files'].append(file_i)
                        continue
                    
                    file_i_path = os.path.join(fs_i[FST._IND_NAME], file_i_title)
                    file_i_id = file_i[FST._IND_ID]
                    
                    if not os.path.exists(file_i_path):
                        local_changes['deleted-files'].append(file_i)
                        continue

                    file_i_mtd = file_i[FST._IND_METADATA]
                    file_i_ModTime_local_old = file_i_mtd['modifiedDateLocal'].replace(tzinfo=tz.tzlocal())
                    file_i_ModTime_local_new = getFileModifiedTime(file_i_path).replace(tzinfo=tz.tzlocal())
                    
                    ##TODO: save and compare file sizes as well
                    
                    if file_i_ModTime_local_old<file_i_ModTime_local_new:
                        local_changes['modified-files'].append(file_i)
                        continue

                fs_children = fs_i[FST._IND_CHILDREN]
                fs_list.extend(fs_children)
                ind_last += 1
                if ind_last>=len(fs_list):
                    break

            fs_tree_new = FST().GetFolderStructure(fs_tree[FST._IND_NAME])
            fs_list = [fs_tree_new]
            ind_last = 0
            while True:
                fs_i = fs_list[ind_last]
                fold_local = fs_i[FST._IND_NAME]
                if fold_local in folders_ignore:
                    local_changes['ignored-folders'].append(fs_i)
                    ind_last += 1
                    if ind_last>=len(fs_list):
                        break
                    continue
                
                if not FST().hasFolder(fs_tree, fold_local):
                    local_changes['new-folders'].append(fs_i)

                    ind_last += 1
                    if ind_last>=len(fs_list):
                        break
                    continue

                fs_files = fs_i[FST._IND_FILES]
                for file_i in fs_files:
                    file_i_title = file_i[FST._IND_NAME]                
                    file_i_path = os.path.join(fs_i[FST._IND_NAME], file_i_title)
                    _file_i_ext_ = os.path.splitext(file_i_title)[1]
                    if _file_i_ext_ in fileTypes_ignore or file_i_title[-1] in fileEndings_ignore:
                        local_changes['ignored-files'].append(file_i)
                        continue
                    if not FST().hasFile(fs_tree, file_i_path):
                        local_changes['new-files'].append(file_i)

                fs_children = fs_i[FST._IND_CHILDREN]
                fs_list.extend(fs_children)
                ind_last += 1
                if ind_last>=len(fs_list):
                    break

        return local_changes_list



    def GetServerChanges(self, dataFileList):
        FST = FolderStructureTree
        server_changes_list = [None]*len(dataFileList)
        for i in range(len(dataFileList)):
            server_changes_list[i] = {'deleted-folders':[], 'deleted-files':[], 
            'modified-folders':[], 'modified-files':[], 
            'new-folders':[], 'new-files':[],
            'renamed-folders':[], 'renamed-files':[],
            'ignored-folders':[], 'ignored-files':[]}
        
        gdfold_str = 'application/vnd.google-apps.folder'
        
        for i_dfl in range(len(dataFileList)):
            start_change_id, folderToSyncDic, fs_tree = dataFileList[i_dfl]
            server_changes = server_changes_list[i_dfl]
            largestChangeId, changes = self.retrieve_all_changes(start_change_id)
            
            folder_sync = folderToSyncDic['foldname']
            gdfolder_id = folderToSyncDic['gdfold-id']
            folders_ignore = folderToSyncDic['folds-ignore']
            fileTypes_ignore = folderToSyncDic['fileTypes-ignore']
            fileEndings_ignore = folderToSyncDic['fileEndings-ignore']

            for ch in changes:
                f_id = None
                if ch['deleted']==True:
                    f_id = ch['fileId']
                else:
                    f_id = ch['file']['id']
                instances = FolderStructureTree().FindByID(fs_tree, f_id)
                if len(instances)>0:
                    #print('Found titles:', [instances[i][FST._IND_NAME] for i in range(len(instances))])
                    assert len(instances)==1
                    f_local = instances[0]
                        
                    if ch['deleted']==True:
                        if f_local[FST._IND_CHILDREN]!=None:
                            server_changes['deleted-folders'].append(f_local)
                        else:
                            server_changes['deleted-files'].append(f_local)
                    
                    else:
                        if ch['file']['mimeType']==gdfold_str:
                            if ch['file']['labels']['trashed']==True:
                                server_changes['deleted-folders'].append(f_local)
                            else:
                                file_i_mtd = f_local[FST._IND_METADATA]
                                file_i_ModTime_server_old = file_i_mtd['modifiedDateServer']
                                file_i_ModTime_server_new = toServerTime(ch['file']['modifiedDate'])
                                
                                if file_i_ModTime_server_old<file_i_ModTime_server_new:
                                    f_title_new = ch['file']['title']
                                    f_modDate = ch['file']['modifiedDate']
                                    server_changes['modified-folders'].append([f_local, f_title_new, f_modDate])
                                
                            
                        else:
                            if ch['file']['labels']['trashed']==True:
                                server_changes['deleted-files'].append(f_local)
                            else:
                                file_i_mtd = f_local[FST._IND_METADATA]
                                file_i_ModTime_server_old = file_i_mtd['modifiedDateServer']
                                file_i_ModTime_server_new = toServerTime(ch['file']['modifiedDate'])
                                
                                if file_i_ModTime_server_old<file_i_ModTime_server_new:
                                    f_title_new = ch['file']['title']
                                    server_changes['modified-files'].append([f_local, f_title_new])

                else: 
                    ##not found
                    if not ch['deleted']:
                        if not ch['file']['labels']['trashed']:
                            parents = ch['file']['parents']
                            for p in parents:
                                p_id = p['id']
                                #print('p_id:', p_id, ch['file']['title'])
                                p_instances = FST().FindByID(fs_tree, p_id)
                                #print([p_instances[i][FST._IND_NAME] for i in range(len(p_instances))])
                                if len(p_instances)>0:
                                    assert len(p_instances)==1
                                    p_local = p_instances[0]
                                    f_title = ch['file']['title']
                                    f_modDate = ch['file']['modifiedDate']
                                    if ch['file']['mimeType']==gdfold_str:
                                        server_changes['new-folders'].append([f_id, f_title, f_modDate, p_local])
                                    else:
                                        server_changes['new-files'].append([f_id, f_title, f_modDate, p_local])
                            
        return server_changes_list                    
            
            
            
            
            


