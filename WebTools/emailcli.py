# emailcli.py


__all__ = ["ImapClient"]

import os
import sys
import imaplib
imaplib._MAXLINE = 200000
import getpass
import email
import datetime
import quopri
from email.message import Message


class ImapClient:
    
    def __init__(self, n_mail=10):
        self.NumberOfMails = n_mail ##number of mails to show
        self.rawMessEncodingType = 'latin-1' ##used for decoding Imap raw message
        self.mailImap = None
        self.debuglevel = 0
        return
        
        
    def setmailserver(self, server):
        self.imapserver = server
        
    def setUserPass(self, user, password):
        self.user = user
        self.password = password
        
    def login(self):
        self.mailImap = imaplib.IMAP4_SSL(self.imapserver)
        self.SetDebug(self.debuglevel)
        try:
            self.mailImap.login(self.user, self.password)
        except imaplib.IMAP4.error:
            print("LOGIN FAILED!!! ")

    def getMailBoxes(self):
        resp, mailboxes = self.mailImap.list()
        if resp == 'OK':
            print("Mailboxes:")
            mailboxes_byline = '\n'.join(map(str, mailboxes))
            print(mailboxes_byline)

    def getMailBoxStatus(self, mailbox_name):
        return self.mailImap.status(mailbox_name, 
            '(MESSAGES RECENT UIDNEXT UIDVALIDITY UNSEEN)')

    def getSelectedMailBoxUnseens_Num(self):
        return self.mailImap.search(None, '(UNSEEN)')

    def getSelectedMailBoxUnseens_UID(self):
        resp, data = self.mailImap.uid('search', None, '(UNSEEN)')
        if resp == 'OK':
            uids = data[0].split()
            return uids
        else:
            return None

    def getSelectedMailBoxDeleted_UID(self):
        resp, data = self.mailImap.uid('search', None, '(Deleted)')
        if resp == 'OK':
            uids = data[0].split()
            return uids
        else:
            return None

    def selectMailbox(self, mailbox):
        resp, data = self.mailImap.select(mailbox)
        if resp == 'OK':
            print("mailbox {} selected".format(mailbox))
            
    def closeMail(self):
        self.mailImap.close()
        
    def logout(self):
        self.mailImap.logout()
        

    def PeekMailbox_All(self, showFullHeader=False, peekbody=False):
        resp, data = self.mailImap.uid('search', None, "ALL")
        if resp != 'OK':
          print("Communication error!")
          return
        
        resp, data_unseen = self.mailImap.uid('search', None, "(UNSEEN)")
        if resp != 'OK':
          print("Communication error!")
          return

        n_read = self.NumberOfMails
        data_split = data[0].split()
        data_unseen_split = data_unseen[0].split()
        n_mail = len(data_split)
        print('number of emails:', len(data_split))
        msg_ids = list(reversed(data_split))[0:min(n_read, n_mail)]
        for m_id in msg_ids:
            unseen = False
            if m_id in data_unseen_split:
                unseen = True
                
            if unseen:
                print('+-'*10)

            m_id = m_id.decode('utf-8')
            print('ID: ', m_id)

            resp, data = self.mailImap.uid('fetch', m_id, '(BODY.PEEK[HEADER])')
            # FLAGS INTERNALDATE ENVELOPE BODY.PEEK[HEADER] BODY.PEEK[TEXT] RFC822 RFC822.SIZE

            if resp != 'OK':
                print("ERROR getting message", m_id)
                return
                                
                    
            msg = email.message_from_string(data[0][1].decode(self.rawMessEncodingType))
            if showFullHeader:
                print(msg)

            print('Raw Date:', msg['Date'])
            date_tuple = email.utils.parsedate_tz(msg['Date'])
            if date_tuple:
                local_date = datetime.datetime.fromtimestamp(email.utils.mktime_tz(date_tuple))
                print("Local Date:", local_date.strftime("%a, %d %b %Y %H:%M:%S"))
                
            subj = self.decodeHeader(msg['Subject'])
            print('Subject : %s' % (subj))
            to__ = self.decodeHeader(msg['To'])
            print('To : %s' % (to__))
            if 'Cc' in msg:
                cc__ = self.decodeHeader(msg['Cc'])
                print('Cc : %s' % (cc__))
            from__ = self.decodeHeader(msg['From'])
            print('From : %s' % (from__))
            if 'Message-ID' in msg:
                message_id__ = self.decodeHeader(msg['Message-ID'])
                print('Message-ID : %s' % (message_id__))
            if 'In-Reply-To' in msg:
                in_reply_to__ = self.decodeHeader(msg['In-Reply-To'])
                print('In-Reply-To : %s' % (in_reply_to__))

            if peekbody:
                resp, data = self.mailImap.uid('fetch', m_id, '(BODY.PEEK[TEXT])')
                if resp != 'OK':
                    print("ERROR getting message", m_id)
                    return
                message_body = data[0][1].decode(self.rawMessEncodingType)
                print("Message PEEK: \n", message_body[0:500])
            
            print("#"*60)
            
    def PeekMailbox_Unseen(self, showFullHeader=False, peekbody=False):
        resp, data = self.mailImap.uid('search', None, "(UNSEEN)")
        if resp != 'OK':
          print("Communication error!")
          return
        
        n_read = self.NumberOfMails
        data_split = data[0].split()
        n_mail = len(data_split)
        print('number of emails:', len(data_split))
        msg_ids = list(reversed(data_split))[0:min(n_read, n_mail)]
        for m_id in msg_ids:
            m_id = m_id.decode('utf-8')
            print('ID: ', m_id)

            resp, data = self.mailImap.uid('fetch', m_id, '(BODY.PEEK[HEADER])')
            # FLAGS INTERNALDATE ENVELOPE BODY.PEEK[HEADER] BODY.PEEK[TEXT] RFC822 RFC822.SIZE

            if resp != 'OK':
                print("ERROR getting message", m_id)
                return
                                
                    
            msg = email.message_from_string(data[0][1].decode(self.rawMessEncodingType))
            if showFullHeader:
                print(msg)

            print('Raw Date:', msg['Date'])
            date_tuple = email.utils.parsedate_tz(msg['Date'])
            if date_tuple:
                local_date = datetime.datetime.fromtimestamp(email.utils.mktime_tz(date_tuple))
                print("Local Date:", local_date.strftime("%a, %d %b %Y %H:%M:%S"))
                
            subj = self.decodeHeader(msg['Subject'])
            print('Subject : %s' % (subj))
            to__ = self.decodeHeader(msg['To'])
            print('To : %s' % (to__))
            if 'Cc' in msg:
                cc__ = self.decodeHeader(msg['Cc'])
                print('Cc : %s' % (cc__))
            from__ = self.decodeHeader(msg['From'])
            print('From : %s' % (from__))
            if 'Message-ID' in msg:
                message_id__ = self.decodeHeader(msg['Message-ID'])
                print('Message-ID : %s' % (message_id__))
            if 'In-Reply-To' in msg:
                in_reply_to__ = self.decodeHeader(msg['In-Reply-To'])
                print('In-Reply-To : %s' % (in_reply_to__))
                
            
            if peekbody:
                resp, data = self.mailImap.uid('fetch', m_id, '(BODY.PEEK[TEXT])')
                if resp != 'OK':
                    print("ERROR getting message", m_id)
                    return
                message_body = data[0][1].decode(self.rawMessEncodingType)
                print("Message PEEK: \n", message_body[0:500])

            print("#"*60)
            
    def decodeHeader(self, header):
        header_decoded_list = email.header.decode_header(header)
        header_dec = ''
        for head_part in header_decoded_list:
            if head_part[1]!=None:
                header_dec = header_dec + head_part[0].decode(head_part[1])
            else:
                #print('head_part[0] :', head_part[0])
                if isinstance(head_part[0], str):
                    header_dec = header_dec + head_part[0]
                else:
                    header_dec = header_dec + head_part[0].decode('utf-8')
        return header_dec
            
            
    def readMailHeader(self, m_id, decodeData=False):
        m_id = m_id.decode('utf-8')

        resp, data = self.mailImap.uid('fetch', m_id, '(BODY.PEEK[HEADER])')
        # FLAGS INTERNALDATE ENVELOPE BODY.PEEK[HEADER] BODY.PEEK[TEXT] RFC822 RFC822.SIZE

        if resp != 'OK':
            print("ERROR getting message", m_id)
            return
                                            
        msg = email.message_from_string(data[0][1].decode(self.rawMessEncodingType))
        if decodeData:
            return None
        else:
            return msg
                
                
    def readMailFlags(self, m_id, decodeData=False):
        m_id = m_id.decode('utf-8')

        resp, data = self.mailImap.uid('fetch', m_id, '(FLAGS)')
        # FLAGS INTERNALDATE ENVELOPE BODY.PEEK[HEADER] BODY.PEEK[TEXT] RFC822 RFC822.SIZE

        if resp != 'OK':
            print("ERROR getting message", m_id)
            return
                                            
        if decodeData:
            return None
        else:
            return data
            
    def readMailEnvelope(self, m_id, decodeData=False):
        m_id = m_id.decode('utf-8')

        resp, data = self.mailImap.uid('fetch', m_id, '(ENVELOPE)')
        # FLAGS INTERNALDATE ENVELOPE BODY.PEEK[HEADER] BODY.PEEK[TEXT] RFC822 RFC822.SIZE

        if resp != 'OK':
            print("ERROR getting message", m_id)
            return
                                            
        if decodeData:
            return None
        else:
            return data
            
    def readMailBody(self, m_id):
        resp, data = self.mailImap.uid('fetch', m_id, '(RFC822)')
        if resp != 'OK':
            print("ERROR getting message", m_id)
            return
        message_body = data[0][1].decode(self.rawMessEncodingType)
        message_body = email.message_from_string(message_body)
        failobj = 'failed'
        char_sets = message_body.get_charsets(failobj=failobj)
        text_parts = []
        html_parts = []
        if message_body.get_content_maintype() == 'multipart':
            ind = -1
            for part in message_body.walk():       
                ind += 1
                if part.get_content_type() == "text/plain":
                    body = part.get_payload(decode=True)
                    if char_sets[ind] != failobj:
                        body = body.decode(char_sets[ind])
                    else:
                        body = body.decode(self.rawMessEncodingType)
                    text_parts.append(body)
                elif part.get_content_type() == "text/html":
                    body = part.get_payload(decode=True)
                    if char_sets[ind] != failobj:
                        body = body.decode(char_sets[ind])
                    else:
                        body = body.decode(self.rawMessEncodingType)
                    html_parts.append(body)
                elif part.get('Content-Disposition')!=None:
                    if part.get('Content-Disposition').startswith("attachment"):
                        fileName = part.get_filename()
                        print("Attachment:", fileName)
                        if bool(fileName):
                            fileName = self.decodeHeader(fileName)
                            detach_dir = '.'
                            if 'attachments' not in os.listdir(detach_dir):
                                os.mkdir('attachments')
                            filePath = os.path.join(detach_dir, 'attachments', fileName)
                            if not os.path.isfile(filePath):
                                print(fileName)
                                fp = open(filePath, 'wb')
                                fp.write(part.get_payload(decode=True))
                                fp.close()                
                else:
                    #part_header = part.items()
                    #print(part_header)
                    pass
        else:
            #print('!!!!!!!!!!!!!! singlepart message !!!!!!!!!!!!')
            #print(char_sets)
            if message_body.get_content_type() == "text/plain":
                body = message_body.get_payload(decode=True)
                if len(char_sets)>0:
                    if char_sets[0] != failobj:
                        body = body.decode(char_sets[0])
                    else:
                        body = body.decode(self.rawMessEncodingType)
                else:
                    body = body.decode(self.rawMessEncodingType)
                text_parts.append(body)
            elif message_body.get_content_type() == "text/html":
                body = message_body.get_payload(decode=True)
                if len(char_sets)>0:
                    if char_sets[0] != failobj:
                        body = body.decode(char_sets[0])
                    else:
                        body = body.decode(self.rawMessEncodingType)
                html_parts.append(body)
            
        return [text_parts, html_parts]
        

    def MarkAsUnread(self, m_id):
        res = self.mailImap.uid('store', m_id, '-FLAGS', r'(\SEEN)')
        return res

    def MarkAsRead(self, m_id):
        res = self.mailImap.uid('store', m_id, '+FLAGS', r'(\SEEN)')
        return res
        
    def DeleteMessage(self, m_id):
        res = self.mailImap.uid('store', m_id, '+FLAGS', r'(\Deleted)')
        return res

    def ExpungeDeletedMessages(self):
        resp, data = self.mailImap.expunge()
        if resp != 'OK':
            print("ERROR trashing messages.")
            return None
        return data
        

    def DeleteMessageUndo(self, m_id):
        res = self.mailImap.uid('store', m_id, '-FLAGS', r'(\Deleted)')
        return res

    def getRawMail(self, m_id):
        resp, data = self.mailImap.uid('fetch', m_id, '(RFC822)')
        if resp != 'OK':
            print("ERROR getting message", m_id)
            return
        return data
        
    def SearchMail(self, charset=None, All=None, Subject=None, 
                    Body=None, From=None):
        search_text = ''
        if All!=None:
            ids_from = self.SearchMail(charset, From=All)
            ids_subj = self.SearchMail(charset, Subject=All)
            ids_body = self.SearchMail(charset, Body=All)
            if ids_from==None or ids_subj==None or ids_body==None:
                return None
            else:
                ids_all = []
                ids_all.extend(ids_from)
                ids_all.extend(ids_subj)
                ids_all.extend(ids_body)
                return sorted(list(set(ids_all)))
        else:
            if From!=None:
                search_text = search_text + '(FROM "%s") '%(From)
            if Subject!=None:
                search_text = search_text + '(SUBJECT "%s") '%(Subject)
            if Body!=None:
                search_text = search_text + '(BODY "%s") '%(Body)
        search_text = search_text.strip()
        #print(search_text)
        resp, data = self.mailImap.uid('search', charset, search_text)
        #print(resp, data)

        if resp != 'OK':
            print("ERROR getting message", m_id)
            return None
        
        data_split = data[0].split()
        return data_split

  

    def SetDebug(self, on=False):
        if on:  
            if self.mailImap:
                self.mailImap.debug = 4
            self.debuglevel = 4
        else:
            if self.mailImap:
                self.mailImap.debug = 0
            self.debuglevel = 0
              
    def IssueUIDCommand(self, command, m_id, args):
        """ Example: command='fetch'
                     args   ='(FLAGS)'
        """
        resp, data = self.mailImap.uid(command, m_id, args)
        return resp, data
  
import os
import sys
import smtplib
# For guessing MIME type based on file name extension
import mimetypes

from optparse import OptionParser

from email import encoders
from email.message import Message
from email.mime.audio import MIMEAudio
from email.mime.base import MIMEBase
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
  
COMMASPACE = ', '
class SMTPEmailer:
    """ Sends SMTP emails
        To: list of receipents (list type or None)
        Cc: list of CCs (list type or None)
    """
    def __init__(self, host, port=587, From=None, To=None, Cc=None, Subject=None, 
            Body=None, In_reply_To=None, Attachments=None):
        """ host: smtp.gmail.com
            port: 587 or 465
        """
        self.host = host
        self.port = port

        self.From = From
        if To!=None:
            self.To = To
        else:
            self.To = []

        self.Subject = Subject
        self.Body = Body

        if Attachments!=None:
            self.Attachments = Attachments
        else:
            self.Attachments = []
            
        if Cc!=None:
            self.Cc = Cc
        else:
            self.Cc = None
        
        if In_reply_To!=None:
            self.In_reply_To = In_reply_To
        else:
            self.In_reply_To = None
            
        self.debuglevel = 0
            
        
    def SetSender(self, sender):
        self.From = sender

    def AddRecipient(self, recepient):
        self.To.append(recepient)
    
    
    def AddCc(self, Cc):
        if self.Cc!=None:
            self.Cc.append(Cc)
        else:
            self.Cc = [Cc]
    
    def SetupMessage(self):    
        msg = MIMEMultipart()
        msg['Subject'] = self.Subject
        msg['To'] = COMMASPACE.join(self.To)
        if self.Cc!=None:
            msg['Cc'] = OMMASPACE.join(self.Cc)
        msg['From'] = self.From
        if self.In_reply_To!=None:
            msg['In-Reply-To'] = self.In_reply_To
        msg.preamble = 'You will not see this in a MIME-aware mail reader.\n'
 
        if self.Body != None:
            # Note: we should handle calculating the charset
            part = MIMEText(self.Body)
            msg.attach(part)
            
        for path in self.Attachments:
            if not os.path.isfile(path):
                print('incorrect attachment: ', path)
                continue
            # Guess the content type based on the file's extension.  Encoding
            # will be ignored, although we should check for simple things like
            # gzip'd or compressed files.
            ctype, encoding = mimetypes.guess_type(path)
            if ctype is None or encoding is not None:
                # No guess could be made, or the file is encoded (compressed), so
                # use a generic bag-of-bits type.
                ctype = 'application/octet-stream'
            maintype, subtype = ctype.split('/', 1)
            part = None
            if maintype == 'text':
                fp = open(path)
                # Note: we should handle calculating the charset
                part = MIMEText(fp.read(), _subtype=subtype)
                fp.close()
            elif maintype == 'image':
                fp = open(path, 'rb')
                part = MIMEImage(fp.read(), _subtype=subtype)
                fp.close()
            elif maintype == 'audio':
                fp = open(path, 'rb')
                part = MIMEAudio(fp.read(), _subtype=subtype)
                fp.close()
            else:
                fp = open(path, 'rb')
                part = MIMEBase(maintype, subtype)
                part.set_payload(fp.read())
                fp.close()
                # Encode the payload using Base64
                encoders.encode_base64(part)
            # Set the filename parameter
            filename = os.path.basename(path)
            part.add_header('Content-Disposition', 'attachment', filename=filename)
            msg.attach(part)
            print(filename, ' attached!')
        # Now send or store the message
        return msg.as_string() 
    
    def Connect(self, username, password):
        server = smtplib.SMTP(self.host, self.port)
        server.set_debuglevel(self.debuglevel)
        server.starttls()
        server.login(username,password)
        return server
    
    def sendEMail(self, username, password):
        try:
            server = self.Connect(username, password)
            message = self.SetupMessage()
            server.sendmail(self.From, self.To, message)
            server.quit()
            return True
        except Exception as e:
            print(str(e))
            return False

    def SetDebug(self, level=0):
        self.debuglevel = level
    
    
    

