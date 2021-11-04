function newPrototypePosition = p_moveTinyBit (prototype)
    moveBy = 0.0000000001;
    newPrototypePosition = prototype + moveBy;
    newPrototypePosition = newPrototypePosition/norm(newPrototypePosition);
end