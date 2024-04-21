for file in sound/*; do
    filename_ext=$(basename "$file")
    filename=$(basename "$file" | cut -f1 -d'.')
    extension=$(echo "$filename_ext" | awk -F. '{print $NF}')
    if [ "$extension" = "wav" ]; then
        ffmpeg -i $file -vn -ar 44100 -ac 2 -b:a 192k sound/"$filename".mp3
        rm $file
    fi
done