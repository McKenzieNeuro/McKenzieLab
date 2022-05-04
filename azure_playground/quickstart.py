"""
https://docs.microsoft.com/en-us/azure/storage/blobs/storage-quickstart-blobs-python?tabs=environment-variable-linux 
"""

import os, uuid, logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

from azure.storage.blob import BlobServiceClient, BlobClient, ContainerClient, __version__


# CONN_STR = "DefaultEndpointsProtocol=https;AccountName=stephenfaystoragetest;AccountKey=WQiLMyK3irl9YkjCJmfD7P5hlU4lNwSk2dimDYYBSKTb57A0DyKgCVHKFfWxsQZEd15O4S2463rU+ASt+Oyo6Q==;EndpointSuffix=core.windows.net"

def create_random_container(conn_str):
    # Create the BlobServiceClient object which will be used to create a container client
    # blob_service_client = BlobServiceClient.from_connection_string(connect_str)
    with BlobServiceClient.from_connection_string(conn_str) as blob_service_client:
        container_name = str(uuid.uuid4()) # random container name
        # create container
        container_client = blob_service_client.create_container(container_name)

def create_container(conn_str:str,container_name:str):
    with BlobServiceClient.from_connection_string(conn_str) as blob_service_client:
        container_client = blob_service_client.create_container(container_name)

def upload_blobs_to_a_container(conn_str:str,container_name:str):
    # Create local directory to hold data
    local_path = "./blob_data"
    if not os.path.exists(local_path): os.mkdir(local_path)

    # Create a file in the local data directory to upload and download
    local_file_name = str(uuid.uuid4()) + ".txt"
    upload_file_path = os.path.join(local_path, local_file_name)

    # Write text to file
    with open(upload_file_path,"w") as file:
        file.write("Hello, World!")

    # Connect to a service client
    with BlobServiceClient.from_connection_string(conn_str) as blob_service_client:
        # Create a blob client using the local file name as the name for the blob
        with blob_service_client.get_blob_client(container=container_name,blob=local_file_name) as blob_client:
            print("Connected initiated blob_client")
            print(f"Uploading to Azure Storage as blob:\n\t{local_file_name}")
            
            # Upload the freshly created file
            with open(upload_file_path,"rb") as data:
                blob_client.upload_blob(data)


# retrieve the connection string stored in environment variable
connect_str  = os.getenv("AZURE_STORAGE_CONNECTION_STRING")
logger.info(f"connect_str={connect_str}")

try:
    print(f"Azure Blob Storave v{__version__} — Python quickstart sample")
    # create_random_container(connect_str)
    upload_blobs_to_a_container(connect_str,"newcontainer")
    
except Exception as ex:
    print("Exception:")
    raise ex


