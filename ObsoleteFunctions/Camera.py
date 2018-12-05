import numpy as np
import cv2

cap = cv2.VideoCapture(1)

print(cap.get(3))
print(cap.get(4))
ret = cap.set(3, cap.get(3)*2)
ret = cap.set(4, cap.get(4)*2)
while(True):
    # Capture frame-by-frame
    ret, frame = cap.read()

    # Our operations on the frame come here
    #gray = cv2.cvtColor(frame, cv2.COLOR_BGR)

    # Display the resulting frame
    cv2.imshow('frame', frame)
    if cv2.waitKey(1) & 0xFF == ord('q'):
        break

# When everything done, release the capture
cap.release()
cv2.destroyAllWindows()